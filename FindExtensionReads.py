import re
import Bio
import os
import sys
import re
import getopt
import pysam
from pysam import AlignmentFile
import Bio
from Bio import SeqIO

class FindExtensionReads(object):
	def __init__(self,roundInput,lastRoundUsedReads,usedReads):
		self.roundInput=roundInput
		self.lastRoundUsedReads=lastRoundUsedReads
		self.usedReads=usedReads
		self.note=''
		self.extensionReads=self.roundInput.elongation.roundDir+"/extensionReads."+self.roundInput.elongation.base.tag+".fa"
		self.log=self.roundInput.elongation.roundDir+"/extensionReads.log"
		if os.path.exists(self.extensionReads)==True and os.path.getsize(self.extensionReads)!=0:
			self.readlog()
		else:
			logfilet=open(self.log,'w')
			self.potentialExtensionReadsAln,self.minimap2Command,self.minimap2Output=self.minimap2()
			logLine='potentialExtensionReadsAln\t'+self.potentialExtensionReadsAln+"\nminimap2Command\t"+self.minimap2Command+"\nminimap2Output\t"+self.minimap2Output+"\n"
			logfilet.writelines(logLine)

			self.minimumExtensionReads()
			logLine='minimumThresholdReadsAln\t'+self.minimumThresholdReadsAln+"\n"
			logLine+='minimumThresholdReadsID\t'+';'.join(self.minimumThresholdReadsID)+"\n"
			logLine+='minimumThresholdExtensionReadsAln\t'+self.minimumThresholdExtensionReadsAln+"\n"
			logLine+='minimumThresholdExtensionReads\t'+self.minimumThresholdExtensionReads+"\n"
			logLine+='minimumThresholdExtensionReadsID\t'+';'.join(self.minimumThresholdExtensionReadsID)+"\n"
			logfilet.writelines(logLine)
			
			
			if self.note=='':
				self.selectMappingQuality=20
				self.selectAlignmentLength=3000
				self.selectNMAlignmentLengthratio=0.1
				
				self.selectReadsNum=0
				self.extensionReadsNum=0
				self.readsExtensionLength=1000
				self.extensionReadsEdge=10
			
				while self.extensionReadsNum<=0 and self.note=='':
					self.selectPotentialExtensionReadsAln=self.roundInput.elongation.roundDir+"/selectPotentialExtensionReads."+self.roundInput.elongation.base.tag+".bam"
					#self.selectPotentialExtensionReadsID=self.samFilter(self.potentialExtensionReadsAln,self.selectPotentialExtensionReadsAln)
					self.selectPotentialExtensionReadsID=self.samFilter(self.minimumThresholdExtensionReadsAln,self.selectPotentialExtensionReadsAln)
					self.selectReadsNum=len(self.selectPotentialExtensionReadsID)

					self.extensionReadsAln=self.roundInput.elongation.roundDir+"/extensionReads."+self.roundInput.elongation.base.tag+".bam"
					self.extensionReads=self.roundInput.elongation.roundDir+"/extensionReads."+self.roundInput.elongation.base.tag+".fa"
					self.extensionReadsID=self.extensionFinder(self.selectPotentialExtensionReadsAln,self.extensionReadsAln,self.extensionReads)
					self.extensionReadsNum=len(self.extensionReadsID)

					if self.extensionReadsNum==0:
						if self.selectMappingQuality==0:
							if self.readsExtensionLength<=10:
								if self.selectAlignmentLength<=500:
									self.extensionReadsEdge=self.extensionReadsEdge+10
									self.selectMappingQuality=20
									self.selectAlignmentLength=3000
									self.readsExtensionLength=1000
								else:
									self.selectAlignmentLength=self.selectAlignmentLength-100
							else:
								self.readsExtensionLength=self.readsExtensionLength-100
						else:
							self.selectMappingQuality=self.selectMappingQuality-10

				logLine='selectPotentialExtensionReadsAln\t'+self.selectPotentialExtensionReadsAln+"\n"
				logLine+='selectPotentialExtensionReadsID\t'+';'.join(self.selectPotentialExtensionReadsID)+"\n"
				logLine+='selectReadsNum\t'+str(self.selectReadsNum)+"\n"
				logLine+='extensionReadsAln\t'+self.extensionReadsAln+"\n"
				logLine+='extensionReads\t'+self.extensionReads+"\n"
				logLine+='extensionReadsID\t'+';'.join(self.extensionReadsID)+"\n"
				logLine+='extensionReadsNum\t'+str(self.extensionReadsNum)+"\n"
				logfilet.writelines(logLine)

			else:
				logLine='selectPotentialExtensionReadsAln\tNone\n'
				logLine+='selectPotentialExtensionReadsID\tNone\n'
				logLine+='selectReadsNum\t0\n'
				logLine+='extensionReadsAln\tNone\n'
				logLine+='extensionReads\tNone\n'
				logLine+='extensionReadsID\tNone\n'
				logLine+='extensionReadsNum\t0\n'
				logfilet.writelines(logLine)
			logLine='selectMappingQuality\t'+str(self.selectMappingQuality)+"\n"
			logLine+='selectAlignmentLength\t'+str(self.selectAlignmentLength)+"\n"
			logLine+='selectNMAlignmentLengthratio\t'+str(self.selectNMAlignmentLengthratio)+"\n"
			logLine+='readsExtensionLength\t'+str(self.readsExtensionLength)+"\n"
			logLine+='extensionReadsEdge\t'+str(self.extensionReadsEdge)+"\n"
			logLine+='note\t'+str(self.note)+"\n"
			logfilet.writelines(logLine)	

			logfilet.close()
	
	def minimumExtensionReads(self):
		self.selectMappingQuality=0
		self.selectAlignmentLength=500
		self.selectNMAlignmentLengthratio=0.1
		self.extensionReadsEdge=self.roundInput.elongation.base.edge
		self.readsExtensionLength=10

		self.minimumThresholdReadsAln=self.roundInput.elongation.roundDir+"/minimumThresholdReads."+self.roundInput.elongation.base.tag+".bam"
		self.minimumThresholdReadsID=self.samFilter(self.potentialExtensionReadsAln,self.minimumThresholdReadsAln)


		self.minimumThresholdExtensionReadsAln=self.roundInput.elongation.roundDir+"/minimumThresholdExtensionReads."+self.roundInput.elongation.base.tag+".bam"
		self.minimumThresholdExtensionReads=self.roundInput.elongation.roundDir+"/minimumThresholdExtensionReads."+self.roundInput.elongation.base.tag+".fa"
		self.minimumThresholdExtensionReadsID=self.extensionFinder(self.minimumThresholdReadsAln,self.minimumThresholdExtensionReadsAln,self.minimumThresholdExtensionReads)

		if len(self.minimumThresholdExtensionReadsID)==0:
			self.note='noExtensionReadsFoundAtMinimumThreshold'
		else:
			self.note=''

	def extensionFinder(self,inputAln,outputAln,outputSeq):
		readslist=[]
		outputSeqFile=open(outputSeq,'w')
		inputAlnFile=AlignmentFile(inputAln,"rb",check_sq=False)
		outputAlnFile=AlignmentFile(outputAln,"wb", template=inputAlnFile)
		#print inputAln,outputAln
		readsExtensionLength=self.readsExtensionLength
		extensionReadsEdge=self.extensionReadsEdge

		for r in inputAlnFile:
			queryID=r.query_name
			query=self.roundInput.elongation.base.readsDict[queryID]
			queryDistance,refDistance,extensionLength,extensionReadsSeq=self.calculateBoundDistance(queryID,query,r)
			if self.readsExtensionLength<=extensionLength:
				if queryDistance<=self.extensionReadsEdge and refDistance<=self.extensionReadsEdge:
					if queryID not in self.usedReads :#or queryID in self.lastRoundUsedReads:
						if queryID not in readslist:
							l='>'+queryID+"\n"+query.seq+"\n"
							outputSeqFile.writelines(l)
							readslist.append(queryID)
					outputAlnFile.write(r)
		
		outputSeqFile.close()
		inputAlnFile.close()
		outputAlnFile.close()
		return readslist
	
	def calculateBoundDistance(self,queryID,query,r):
		if r.is_reverse:
			queryseq=query.seq.reverse_complement()
		else:
			queryseq=query.seq
		
		qs,qe,rs,re=self.findAlnPosition(queryseq,self.roundInput.inputSeedSequence,r)
		if self.roundInput.elongation.base.flag=='left':
			queryDistance=qs
			extensionReadsSeq=queryseq
			refDistance=len(self.roundInput.inputSeedSequence.seq)-re
			extensionLength=len(queryseq)-qe
		else:
			queryDistance=len(queryseq)-qe
			extensionReadsSeq=queryseq
			refDistance=rs
			extensionLength=qs
		if qs<0 or len(queryseq)<qe or qe<0 or rs<0 or re<0 :#or len(self.roundInput.inputSeedSequence.seq)<re:
			print ('wrong seqend',qs,qe,rs,re)
			print (queryDistance,refDistance,extensionLength,extensionReadsSeq)
			print (qs<0)
			print (len(queryseq)<qe)
			print (qe<0)
			print (rs<0)
			print (re<0)
			print (len(self.roundInput.inputSeedSequence.seq)<re)
			print (len(self.roundInput.inputSeedSequence.seq),re)
			sys.exit()
		return queryDistance,refDistance,extensionLength,extensionReadsSeq
	
	def findAlnPosition(self,queryseq,refseq,r):
		ct=r.cigartuples
		scpstart=ct[0]
		if scpstart[0]==4:
			scps=scpstart[1]
		else:
			scps=0
	
		scpend=ct[-1]
		if scpend[0]==4:
			scpe=scpend[1]
		else:
			scpe=0

		if queryseq==r.query_sequence or r.query_alignment_sequence==None:
			qs=r.query_alignment_start
			qe=r.query_alignment_end
		else:
			qaln=r.query_alignment_sequence
			qs=str(queryseq).index(qaln)
			qe=qs+len(r.query_alignment_sequence)
		rs=r.reference_start
		re=r.reference_end
		if qs<0:
			print ('wrong qs',qs)
		if qe>len(queryseq):
			print ('wrong qe',qe,len(queryseq))
		return qs,qe,rs,re

		

	def samFilter(self,inputAln,outputAln):
		samFile=AlignmentFile(inputAln,"rb",check_sq=False)
		outputBamFile=AlignmentFile(outputAln,"wb", template=samFile)
		
		readslist=[]
		print (inputAln)
		sMQ=self.selectMappingQuality
		sAlignmentLength=self.selectAlignmentLength
		sNMAlignmentLengthr=self.selectNMAlignmentLengthratio

		for r in samFile:
			if r.is_unmapped==False and r.mapping_quality>=sMQ:
				for i in r.tags:
					if i[0]=='NM':
						NM=i[1]
				AlignmentLength=len(r.query_alignment_sequence)
				if  AlignmentLength>=sAlignmentLength and float(NM)/AlignmentLength<=sNMAlignmentLengthr:
					outputBamFile.write(r)
					if r.query_name not in readslist:
						readslist.append(r.query_name)
		samFile.close()
		outputBamFile.close()
		return readslist
	
	def jellyfish(self, kmerseq, kmer_size, kmer_num, min_count=1):
		"""
		调用jellyfish命令，统计kmerseq的kmer_size长度的kmer的出现次数
		按照kmer的出现次数，从低到高排序
		从上到下，取kmer_num个低频kmer（但忽略出现次数低于min_count的kmer）
		每行为一个kmer
		将kmer_num个kmer保存为self.roundInput.elongation.roundDir + '/kmers.txt'
		
		参数:
			kmerseq: 需要统计kmer的序列
			kmer_size: kmer的长度
			kmer_num: 需要提取的kmer数量
			min_count: 最小出现次数阈值，低于此值的kmer将被忽略（默认为1）
		
		返回:
			保存kmer的文件路径
		"""
		import subprocess
		import os
		import tempfile
		
		# 创建输出目录（如果不存在）
		output_dir = self.roundInput.elongation.roundDir
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
		
		# 按count升序排序（从低到高）
		kmers.sort(key=lambda x: x[1])
		
		# 取前kmer_num个kmer写入输出文件
		with open(kmer_output, 'w') as f:
			for i, (kmer, count) in enumerate(kmers):
				if i >= kmer_num:
					break
				f.write(f"{kmer}\n")
				print(f"低频Kmer {i+1}: {kmer} (count: {count})")
		
		# 清理临时文件
		os.remove(temp_fasta)
		os.remove(f"{output_dir}/all_kmers.txt")
		if os.path.exists(jf_db):
			os.remove(jf_db)
		
		return str(kmer_output)

	def seqkit(self, kmer_list, reads):
		"""
		调用seqkit grep命令从reads中提取包含kmer_list中指定k-mer的序列
		
		参数:
			kmer_list: 包含k-mer的文件路径
			reads: 需要搜索的reads文件路径
		
		返回:
			输出文件的路径
		"""
		import subprocess
		import os
		
		# 创建输出目录（如果不存在）
		output_dir = self.roundInput.elongation.roundDir
		if not os.path.exists(output_dir):
			os.makedirs(output_dir)
		
		# 设置输出文件路径
		output_file = os.path.join(output_dir, "seqkitOutput.fa")
		
		# 执行seqkit grep命令
		# -s: 匹配序列，而不是ID
		# -m: 模式匹配
		# -f: 指定模式文件
		cmd = f"seqkit grep  -f {kmer_list}  -s {reads} > {output_file}"
		
		# 执行命令
		process = subprocess.run(cmd, shell=True, check=True)
		
		# 检查命令是否成功执行
		if process.returncode != 0:
			raise Exception(f"seqkit命令执行失败，返回码: {process.returncode}")
		
		return str(output_file)

	def blast(self, kmer_list, reads):
		"""
		使用blastn命令比对kmer_list与reads，并提取匹配的reads序列
		
		参数:
			kmer_list: 包含k-mer的文件路径
			reads: 需要搜索的reads文件路径
		
		返回:
			输出文件的路径
		"""
		import subprocess
		import os
		from Bio import SeqIO
		
		# 创建输出目录（如果不存在）
		output_dir = self.roundInput.elongation.roundDir
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
		
		# 执行blastn命令
		cmd = f"blastn -query {kmer_fasta} -subject {reads} -outfmt 6 -out {blast_out}"
		process = subprocess.run(cmd, shell=True, check=True)
		
		# 检查命令是否成功执行
		if process.returncode != 0:
			raise Exception(f"blastn命令执行失败，返回码: {process.returncode}")
		
		# 从blast结果中提取匹配的reads ID
		matched_ids = set()
		with open(blast_out, 'r') as f:
			for line in f:
				fields = line.strip().split('\t')
				if len(fields) >= 2:
					matched_ids.add(fields[1])  # 第二列是subject ID
		
		# 从reads文件中提取匹配的序列
		with open(output_file, 'w') as out_f:
			for record in SeqIO.parse(reads, "fasta"):
				if record.id in matched_ids:
					SeqIO.write(record, out_f, "fasta")
		
		# 清理临时文件
		os.remove(kmer_fasta)
		
		return str(output_file)

	def kmerfilter(self, inputSeq, length, kmer_size, kmer_num):
		"""
		使用k-mer过滤方法筛选HiFi读段
		
		参数:
			inputSeq: 输入序列文件
			length: 长度阈值
			kmer_size: k-mer的大小
			kmer_num: 需要的k-mer数量
		
		返回:
			过滤后的reads文件路径
		"""
		import os
		import subprocess
		import math
		
		# 获取平均读段长度
		kmerseqlength = 0
		hifi_stat = self.roundInput.elongation.base.out + "/HiFi.reads.stat"
		
		if os.path.exists(hifi_stat) and os.path.getsize(hifi_stat) > 0:
			# 从HiFi.reads.stat文件读取数据
			number = 0
			total_len = 0
			with open(hifi_stat, 'r') as f:
				for line in f:
					parts = line.strip().split('\t')
					if parts[0] == "Number":
						number = int(parts[1])
					elif parts[0] == "TolalLenth":
						total_len = int(parts[1])
		
			if number > 0:
				kmerseqlength = total_len / number
		else:
			# 使用awk统计读段平均长度
			reads_file = self.roundInput.elongation.base.reads
			cmd = f"awk '{{if(NR%4==2) sum+=length($0); count++}} END {{print sum/count}}' {reads_file}"
			if reads_file.endswith('.fasta') or reads_file.endswith('.fa'):
				cmd = f"awk '{{if(NR%2==0) sum+=length($0); count++}} END {{print sum/count}}' {reads_file}"
			
			result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
			if result.returncode == 0 and result.stdout.strip():
				kmerseqlength = float(result.stdout.strip())
		
		# 获取DEGAP.py的参数
		is_right_flag = self.roundInput.elongation.base.flag == 'right'
		
		# 根据flag参数提取序列
		kmerseq = ""
		for seq_record in SeqIO.parse(self.roundInput.inputSeq, "fasta"):
			seq = str(seq_record.seq)
			target_len = int(kmerseqlength * 0.1)  # 取10%的长度
			
			if is_right_flag:
				# 从左端提取序列
				kmerseq = seq[:min(target_len, len(seq))]
			else:
				# 从右端提取序列
				kmerseq = seq[max(0, len(seq) - target_len):]
			
			# 只处理第一条序列
			break
		
		# 使用jellyfish生成k-mer列表
		kmer_list = self.jellyfish(kmerseq, kmer_size, kmer_num)
		
		# 使用seqkit筛选包含这些k-mer的reads
		seqkitOutput = self.seqkit(kmer_list, self.roundInput.elongation.base.reads)
		
		# 使用blast进一步过滤
		blastOutput = self.blast(kmer_list, seqkitOutput)
		if os.path.exists(blastOutput) and os.path.getsize(blastOutput) > 0:
			return blastOutput
		else:
			return seqkitOutput

	def minimap2(self):
		alnname=self.roundInput.elongation.roundDir+"/potentialExtensionReads."+self.roundInput.elongation.base.tag+".bam"
		alnname1=self.roundInput.elongation.roundDir+"/extensionReads."+self.roundInput.elongation.base.tag+".fa"
		
		# 如果目标文件已存在，先删除它以避免追加写入
		if os.path.exists(alnname):
			try:
				os.remove(alnname)
			except Exception as e:
				print(f"警告：无法删除已存在的文件 {alnname}: {e}")
		
		#新增滤网层
		if os.path.getsize(self.roundInput.elongation.base.reads) > 1024*1024*1024:  # 如果文件大于1GB
			print("输入文件大于1GB，应用k-mer过滤...")
			filtered_reads = self.kmerfilter(self.roundInput.inputSeq, length=100, kmer_size=51, kmer_num=10)
			# 确保过滤结果有效
			if os.path.exists(filtered_reads) and os.path.getsize(filtered_reads) > 0:
				reads_to_use = filtered_reads
				print(f"使用过滤后的reads文件: {reads_to_use}")
			else:
				reads_to_use = self.roundInput.elongation.base.reads
				print(f"k-mer过滤未生成有效文件，使用原始reads: {reads_to_use}")
		else:
			reads_to_use = self.roundInput.elongation.base.reads

		commandline="minimap2 -t "+self.roundInput.elongation.base.thread+" -Y -ax asm20 "+self.roundInput.inputSeq+" "+reads_to_use+" | samtools view -bS >"+alnname
		
		# 如果文件已存在且有效，直接返回
		if os.path.exists(alnname)==True and os.path.getsize(alnname)!=0 and os.path.exists(alnname1)==True and os.path.getsize(alnname1)!=0:
			return alnname,commandline,str(0)
		else:
			import subprocess
			import time
			
			# 使用subprocess库替代os.system，增加超时控制
			try:
				print(f"执行命令: {commandline}")
				start_time = time.time()
				# 设置超时时间为30分钟
				result = subprocess.run(commandline, shell=True, timeout=1800, capture_output=True)
				
				# 检查返回代码
				if result.returncode != 0:
					print(f"minimap2命令执行失败，返回码: {result.returncode}")
					print(f"错误输出: {result.stderr.decode('utf-8', errors='ignore')}")
					
					# 重试最多2次
					retries = 0
					while result.returncode != 0 and retries < 2:
						retries += 1
						print(f"第{retries}次重试minimap2命令...")
						result = subprocess.run(commandline, shell=True, timeout=1800, capture_output=True)
						
						if result.returncode == 0:
							break
							
						if retries == 2:
							print("minimap2多次重试失败，无法正确执行!")
							sys.exit(1)
				
				# 检查输出文件大小
				if os.path.exists(alnname):
					file_size = os.path.getsize(alnname)
					if file_size == 0:
						print(f"警告：minimap2生成的BAM文件大小为0字节")
					else:
						print(f"minimap2生成的BAM文件大小: {file_size} 字节")
						
				return alnname, commandline, str(result.returncode)
				
			except subprocess.TimeoutExpired:
				print("minimap2命令执行超时，可能是处理大文件导致")
				# 如果超时，尝试终止进程并清理可能的部分文件
				if os.path.exists(alnname):
					os.remove(alnname)
				sys.exit(1)
			except Exception as e:
				print(f"执行minimap2命令时发生错误: {e}")
				sys.exit(1)

	def readlog(self):
		logfilet=open(self.log,'r')
		for row in logfilet:
			row1=row.rstrip().split('\t')
			if row1[0]=='potentialExtensionReadsAln':
				self.potentialExtensionReadsAln=row1[1]
			elif row1[0]=='minimap2Command':
				self.minimap2Command=row1[1]
			elif row1[0]=='minimap2Output':
				self.minimap2Output=row1[1]
			elif row1[0]=='minimumThresholdReadsAln':
				self.minimumThresholdReadsAln=row1[1]
			elif row1[0]=='minimumThresholdReadsID':
				self.minimumThresholdReadsID=row1[1].split(';')
			elif row1[0]=='minimumThresholdExtensionReadsAln':
				self.minimumThresholdExtensionReadsAln=row1[1]
			elif row1[0]=='minimumThresholdExtensionReads':
				self.minimumThresholdExtensionReads=row1[1]
			elif row1[0]=='minimumThresholdExtensionReadsID':
				self.minimumThresholdExtensionReadsID=row1[1].split(';')
			elif row1[0]=='selectPotentialExtensionReadsAln':
				self.selectPotentialExtensionReadsAln=row1[1]

			elif row1[0]=='selectPotentialExtensionReadsID':
				if row1[1]!='None':
					self.selectPotentialExtensionReadsID=row1[1].split(';')
				else:
					self.selectPotentialExtensionReadsID=[]
			elif row1[0]=='selectReadsNum':
				self.selectReadsNum=int(row1[1])
			elif row1[0]=='extensionReadsAln':
				self.extensionReadsAln=row1[1]
			elif row1[0]=='extensionReads':
				self.extensionReads=row1[1]
			elif row1[0]=='extensionReadsID':
				if row1[1]!='None':
					self.extensionReadsID=row1[1].split(';')
				else:
					self.extensionReadsID=[]
			elif row1[0]=='extensionReadsNum':
				self.extensionReadsNum=int(row1[1])
			elif row1[0]=='selectMappingQuality':
				self.selectMappingQuality=int(row1[1])
			elif row1[0]=='selectAlignmentLength':
				self.selectAlignmentLength=int(row1[1])
			elif row1[0]=='selectNMAlignmentLengthratio':
				self.selectNMAlignmentLengthratio=float(row1[1])
			elif row1[0]=='readsExtensionLength':
				self.readsExtensionLength=int(row1[1])
			elif row1[0]=='extensionReadsEdge':
				self.extensionReadsEdge=int(row1[1])
			elif row1[0]=='note':
				if len(row1)==1:
					self.note=''
				else:
					self.note=row1[1]
		logfilet.close()
