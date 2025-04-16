#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import sys
import subprocess
import argparse
from Bio import SeqIO
import math
import tempfile
import time
import shutil

def jellyfish(kmerseq, kmer_size, kmer_num, output_dir, min_count=2):
    """
    调用jellyfish命令，统计kmerseq的kmer_size长度的kmer的出现次数
    按照kmer的出现次数，从低到高排序，取kmer_num个低频但非唯一的k-mer
    
    参数:
        kmerseq: 需要统计kmer的序列
        kmer_size: kmer的长度
        kmer_num: 需要提取的k-mer数量
        output_dir: 输出目录
        min_count: 最小k-mer计数，低于此值的k-mer将被忽略(默认为2，可避免测序错误)
    
    返回:
        保存kmer的文件路径
    """
    # 创建输出目录（如果不存在）
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    # 创建临时文件保存输入序列
    temp_fasta = os.path.join(output_dir, "temp_input.fa")
    with open(temp_fasta, 'w') as f:
        f.write(">temp_seq\n")
        f.write(kmerseq)
    
    # 设置输出文件路径
    jf_db = os.path.join(output_dir, "kmer_counts.jf")
    kmer_output = os.path.join(output_dir, "kmers.txt")
    
    # 调用jellyfish count统计kmer
    print(f"[INFO] 运行jellyfish count统计{kmer_size}-mers...")
    count_cmd = f"jellyfish count -m {kmer_size} -s 100M -o {jf_db} {temp_fasta}"
    subprocess.run(count_cmd, shell=True, check=True)
    
    # 使用jellyfish dump导出所有kmer及其计数
    dump_cmd = f"jellyfish dump -c {jf_db} > {output_dir}/all_kmers.txt"
    subprocess.run(dump_cmd, shell=True, check=True)
    
    # 读取所有kmer，按计数排序并获取前kmer_num个
    kmers = []
    with open(f"{output_dir}/all_kmers.txt", 'r') as f:
        for line in f:
            parts = line.strip().split()
            if len(parts) == 2:
                kmer, count = parts
                count = int(count)
                # 仅保留计数大于或等于min_count的k-mer
                if count >= min_count:
                    kmers.append((kmer, count))
    
    # 按count升序排序(从低到高)
    kmers.sort(key=lambda x: x[1])
    
    # 取前kmer_num个低频k-mer写入输出文件
    with open(kmer_output, 'w') as f:
        for i, (kmer, count) in enumerate(kmers):
            if i >= kmer_num:
                break
            f.write(f"{kmer}\n")
            print(f"[INFO] 低频Kmer {i+1}: {kmer} (count: {count})")
    
    # 添加统计信息
    if len(kmers) < kmer_num:
        print(f"[WARNING] 只找到{len(kmers)}个满足条件的k-mer，少于请求的{kmer_num}个")
    
    # 清理临时文件
    os.remove(temp_fasta)
    os.remove(f"{output_dir}/all_kmers.txt")
    if os.path.exists(jf_db):
        os.remove(jf_db)
    
    print(f"[INFO] 已生成包含{min(kmer_num, len(kmers))}个低频k-mers的文件: {kmer_output}")
    return kmer_output

def seqkit(kmer_list, reads, output_dir):
    """
    调用seqkit grep命令从reads中提取包含kmer_list中指定k-mer的序列
    
    参数:
        kmer_list: 包含k-mer的文件路径
        reads: 需要搜索的reads文件路径
        output_dir: 输出目录
    
    返回:
        输出文件的路径
    """
    # 创建输出目录（如果不存在）
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    # 设置输出文件路径
    output_file = os.path.join(output_dir, "seqkitOutput.fa")
    
    # 执行seqkit grep命令
    # -f: 指定模式文件
    # -s: 匹配序列，而不是ID
    print(f"[INFO] 运行seqkit grep提取包含k-mer的序列...")
    cmd = f"seqkit grep -f {kmer_list} -s {reads} > {output_file}"
    
    try:
        # 记录开始时间
        start_time = time.time()
        # 执行命令
        process = subprocess.run(cmd, shell=True, check=True)
        end_time = time.time()
        
        # 检查输出文件
        if os.path.exists(output_file) and os.path.getsize(output_file) > 0:
            # 计算序列数量
            seq_count = 0
            for _ in SeqIO.parse(output_file, "fasta"):
                seq_count += 1
            
            print(f"[INFO] seqkit提取成功，得到{seq_count}条序列，用时{end_time-start_time:.2f}秒")
            print(f"[INFO] 输出文件: {output_file}")
            return output_file
        else:
            print(f"[WARNING] seqkit未能提取到任何序列")
            return None
    except Exception as e:
        print(f"[ERROR] seqkit命令执行失败: {e}")
        return None

def blast(kmer_list, reads, output_dir):
    """
    使用blastn命令比对kmer_list与reads，并提取匹配的reads序列
    
    参数:
        kmer_list: 包含k-mer的文件路径
        reads: 需要搜索的reads文件路径
        output_dir: 输出目录
    
    返回:
        输出文件的路径
    """
    # 创建输出目录（如果不存在）
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    # 设置输出文件路径
    blast_out = os.path.join(output_dir, "blastOutput.txt")
    output_file = os.path.join(output_dir, "blastOutput.fa")
    
    # 创建包含kmer的FASTA文件（blastn需要FASTA格式输入）
    kmer_fasta = os.path.join(output_dir, "kmers.fasta")
    with open(kmer_list, 'r') as kmer_in, open(kmer_fasta, 'w') as kmer_out:
        kmer_count = 0
        for line in kmer_in:
            kmer = line.strip()
            if kmer:
                kmer_out.write(f">kmer_{kmer_count}\n{kmer}\n")
                kmer_count += 1
    
    # 确认至少有一个k-mer
    if kmer_count == 0:
        print("[ERROR] k-mer列表为空，无法执行BLAST")
        return None
    
    # 执行blastn命令
    print(f"[INFO] 运行blastn比对{kmer_count}个k-mers...")
    cmd = f"blastn -query {kmer_fasta} -subject {reads} -outfmt 6 -out {blast_out}"
    
    try:
        # 记录开始时间
        start_time = time.time()
        # 执行命令
        process = subprocess.run(cmd, shell=True, check=True)
        end_time = time.time()
        
        # 检查命令是否成功执行
        if process.returncode != 0:
            print(f"[ERROR] blastn命令执行失败，返回码: {process.returncode}")
            return None
        
        # 从blast结果中提取匹配的reads ID
        matched_ids = set()
        with open(blast_out, 'r') as f:
            for line in f:
                fields = line.strip().split('\t')
                if len(fields) >= 2:
                    matched_ids.add(fields[1])  # 第二列是subject ID
        
        print(f"[INFO] BLAST找到{len(matched_ids)}个匹配的序列，用时{end_time-start_time:.2f}秒")
        
        # 如果没有匹配结果，返回None
        if len(matched_ids) == 0:
            print("[WARNING] 未找到任何匹配的序列")
            return None
        
        # 从reads文件中提取匹配的序列
        count = 0
        with open(output_file, 'w') as out_f:
            for record in SeqIO.parse(reads, "fasta"):
                if record.id in matched_ids:
                    SeqIO.write(record, out_f, "fasta")
                    count += 1
        
        print(f"[INFO] 已将{count}条匹配序列写入输出文件: {output_file}")
        
        # 清理临时文件
        os.remove(kmer_fasta)
        
        return output_file
    
    except Exception as e:
        print(f"[ERROR] 执行BLAST时发生错误: {e}")
        return None

def kmerfilter(seqright=None, seqleft=None, reads=None, output_dir="./output", flag="left", kmer_size=31, kmer_num=10, extract_length=1000):
    """
    使用k-mer过滤方法筛选HiFi读段
    
    参数:
        seqright: 右端序列文件路径
        seqleft: 左端序列文件路径
        reads: 原始reads文件
        output_dir: 输出目录
        flag: 选择从哪端提取序列 ('left'或'right')
        kmer_size: k-mer的大小
        kmer_num: 需要的k-mer数量
        extract_length: 从序列末端提取的长度
    
    返回:
        过滤后的reads文件路径
    """
    # 确保输出目录存在
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    print("[INFO] 开始k-mer过滤流程...")
    print(f"[INFO] 参数：flag={flag}, kmer_size={kmer_size}, kmer_num={kmer_num}")
    
    # 根据flag参数确定使用哪个序列文件并提取相应区域
    input_seq = None
    kmerseq = ""
    
    if flag == "left":
        if not seqleft:
            print("[ERROR] flag=left但未提供--seqleft参数")
            return None
        print(f"[INFO] 使用seqleft文件: {seqleft}")
        try:
            # 从左端文件的右侧(3'端)提取序列
            for seq_record in SeqIO.parse(seqleft, "fasta"):
                seq = str(seq_record.seq)
                if len(seq) <= extract_length:
                    kmerseq = seq
                    print(f"[INFO] seqleft序列长度小于{extract_length}bp，使用全长: {len(seq)}bp")
                else:
                    kmerseq = seq[-extract_length:]
                    print(f"[INFO] 从seqleft右端(3'端)提取{extract_length}bp序列用于k-mer生成")
                break  # 只处理第一条序列
        except Exception as e:
            print(f"[ERROR] 读取seqleft文件时出错: {e}")
            return None
    elif flag == "right":
        if not seqright:
            print("[ERROR] flag=right但未提供--seqright参数")
            return None
        print(f"[INFO] 使用seqright文件: {seqright}")
        try:
            # 从右端文件的左侧(5'端)提取序列
            for seq_record in SeqIO.parse(seqright, "fasta"):
                seq = str(seq_record.seq)
                if len(seq) <= extract_length:
                    kmerseq = seq
                    print(f"[INFO] seqright序列长度小于{extract_length}bp，使用全长: {len(seq)}bp")
                else:
                    kmerseq = seq[:extract_length]
                    print(f"[INFO] 从seqright左端(5'端)提取{extract_length}bp序列用于k-mer生成")
                break  # 只处理第一条序列
        except Exception as e:
            print(f"[ERROR] 读取seqright文件时出错: {e}")
            return None
    else:
        print(f"[ERROR] 不支持的flag值: {flag}，请使用'left'或'right'")
        return None
    
    if not kmerseq:
        print("[ERROR] 无法提取序列")
        return None
    
    # 保存提取的序列以便查看
    extracted_seq_file = os.path.join(output_dir, "extracted_sequence.fa")
    with open(extracted_seq_file, 'w') as f:
        f.write(f">extracted_{flag}_end\n{kmerseq}\n")
    print(f"[INFO] 提取的序列已保存到: {extracted_seq_file}")
    
    # 使用jellyfish生成k-mer列表
    kmer_list = jellyfish(kmerseq, kmer_size, kmer_num, output_dir)
    
    # 使用seqkit筛选包含这些k-mer的reads
    seqkitOutput = seqkit(kmer_list, reads, output_dir)
    
    # 如果seqkit失败，返回None
    if not seqkitOutput:
        print("[ERROR] seqkit筛选失败")
        return None
    
    # 使用blast进一步过滤
    blastOutput = blast(kmer_list, seqkitOutput, output_dir)
    
    # 检查最终输出
    final_output = os.path.join(output_dir, "filter.fa")
    if blastOutput and os.path.exists(blastOutput) and os.path.getsize(blastOutput) > 0:
        print(f"[INFO] k-mer过滤完成，最终输出文件: {blastOutput}")
        # 将blast结果复制为统一的输出文件
        shutil.copy2(blastOutput, final_output)
        print(f"[INFO] 最终过滤结果已保存至: {final_output}")
        return final_output
    else:
        print(f"[WARNING] BLAST过滤失败或结果为空，使用seqkit的输出: {seqkitOutput}")
        # 将seqkit结果复制为统一的输出文件
        shutil.copy2(seqkitOutput, final_output)
        print(f"[INFO] 最终过滤结果已保存至: {final_output}")
        return final_output

def main():
    """
    主函数，处理命令行参数并执行相应的功能
    """
    parser = argparse.ArgumentParser(description="K-mer过滤工具集，包含jellyfish提取k-mer、seqkit筛选和BLAST比对功能")
    
    subparsers = parser.add_subparsers(dest="command", help="要执行的命令")
    
    # jellyfish子命令
    jellyfish_parser = subparsers.add_parser("jellyfish", help="使用jellyfish提取k-mer")
    jellyfish_parser.add_argument("-i", "--input", required=True, help="输入序列(字符串或文件路径)")
    jellyfish_parser.add_argument("-k", "--kmer_size", type=int, default=31, help="k-mer大小")
    jellyfish_parser.add_argument("-n", "--kmer_num", type=int, default=10, help="提取的k-mer数量")
    jellyfish_parser.add_argument("-o", "--output_dir", default="./output", help="输出目录")
    jellyfish_parser.add_argument("-f", "--from_file", action="store_true", help="从文件读取输入序列而不是直接使用字符串")
    
    # seqkit子命令
    seqkit_parser = subparsers.add_parser("seqkit", help="使用seqkit筛选包含指定k-mer的reads")
    seqkit_parser.add_argument("-k", "--kmer_list", required=True, help="包含k-mer的文件路径")
    seqkit_parser.add_argument("-r", "--reads", required=True, help="需要搜索的reads文件路径")
    seqkit_parser.add_argument("-o", "--output_dir", default="./output", help="输出目录")
    
    # blast子命令
    blast_parser = subparsers.add_parser("blast", help="使用BLAST比对k-mer和reads")
    blast_parser.add_argument("-k", "--kmer_list", required=True, help="包含k-mer的文件路径")
    blast_parser.add_argument("-r", "--reads", required=True, help="需要搜索的reads文件路径")
    blast_parser.add_argument("-o", "--output_dir", default="./output", help="输出目录")
    
    # kmerfilter子命令 - 新的参数格式
    kmerfilter_parser = subparsers.add_parser("kmerfilter", help="综合使用k-mer方法过滤reads")
    kmerfilter_parser.add_argument("--seqright", help="右端序列文件路径")
    kmerfilter_parser.add_argument("--seqleft", help="左端序列文件路径")
    kmerfilter_parser.add_argument("--reads", required=True, help="需要过滤的reads文件路径")
    kmerfilter_parser.add_argument("--flag", choices=["left", "right"], default="left", help="选择从哪端提取序列(left或right)")
    kmerfilter_parser.add_argument("-k", "--kmer_size", type=int, default=31, help="k-mer大小")
    kmerfilter_parser.add_argument("-n", "--kmer_num", type=int, default=10, help="需要的k-mer数量")
    kmerfilter_parser.add_argument("-l", "--extract_length", type=int, default=1000, help="从序列末端提取的长度")
    kmerfilter_parser.add_argument("-o", "--output_dir", default="./output", help="输出目录")
    
    args = parser.parse_args()
    
    # 根据命令执行相应的功能
    if args.command == "jellyfish":
        if args.from_file:
            # 从文件读取序列
            try:
                with open(args.input, 'r') as f:
                    kmerseq = ""
                    for record in SeqIO.parse(f, "fasta"):
                        kmerseq = str(record.seq)
                        break
                if not kmerseq:
                    print("[ERROR] 无法从输入文件中读取序列")
                    return 1
            except Exception as e:
                print(f"[ERROR] 读取输入文件时出错: {e}")
                return 1
        else:
            # 直接使用输入字符串
            kmerseq = args.input
        
        jellyfish(kmerseq, args.kmer_size, args.kmer_num, args.output_dir)
    
    elif args.command == "seqkit":
        seqkit(args.kmer_list, args.reads, args.output_dir)
    
    elif args.command == "blast":
        blast(args.kmer_list, args.reads, args.output_dir)
    
    elif args.command == "kmerfilter":
        kmerfilter(
            seqright=args.seqright,
            seqleft=args.seqleft,
            reads=args.reads,
            output_dir=args.output_dir,
            flag=args.flag,
            kmer_size=args.kmer_size,
            kmer_num=args.kmer_num,
            extract_length=args.extract_length
        )
    
    else:
        parser.print_help()
        return 1
    
    return 0

if __name__ == "__main__":
    sys.exit(main()) 