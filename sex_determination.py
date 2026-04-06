# coding: utf-8

# Author: Julien Buratti
#!/usr/bin/env python3

import os
import glob
import pysam
import sys

sample_l = []
sex_d = {}
count_d = {} # 1. SRY count를 저장할 딕셔너리 추가

bams = glob.glob("*/*.dedup.bam")

print("\n************************************")
print("Sex determination script opening.")
print("************************************\n")

for i in bams:
    sample = os.path.basename(i).split(".")[0]
    sample_l.append(sample)
    bamfile = pysam.AlignmentFile(i, "rb", check_sq=False)
    sry_count = bamfile.count(contig='chrY', start=2786989, stop=2787603, until_eof=False, read_callback='all')
    
    print(f"[{sample}] SRY count: {sry_count}")

    if sry_count >= 50:
        sex = "M"
    elif sry_count <= 10:
        sex = "F"
    else:
        sex = "?"
    
    sex_d[sample] = sex
    count_d[sample] = sry_count # 2. 계산된 SRY count를 딕셔너리에 저장

# 3. samples.txt 생성 블록
if not os.path.isfile("samples.txt"):
    with open('samples.txt', 'w') as samp:
        samp.write("sample\tsex\n")
        for s in sample_l:
            samp.write(s + "\t" + sex_d[s] + "\n")
    print("samples.txt generated.")
else:
    print("samples.txt already exists. No changes.")

# 4. sry_counts.txt 생성 블록 (추가된 부분)
if not os.path.isfile("sry_counts.txt"):
    with open('sry_counts.txt', 'w') as sry_out:
        # 헤더 작성: sample, sex, sry_count
        sry_out.write("sample\tsex\tsry_count\n")
        
        for s in sample_l:
            # 기존 정보와 함께 count_d[s]에 저장해둔 SRY count 값을 기록
            sry_out.write(f"{s}\t{sex_d[s]}\t{count_d[s]}\n")
            
    print("sry_counts.txt generated.")
else:
    print("sry_counts.txt already exists. No changes.")

print(f"Sex determination done for {len(sample_l)} samples.")
print("\nSex determination script job done!\n")
