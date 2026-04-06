# coding: utf-8

# Author: Julien Buratti
#!/usr/bin/env python3

import os
import glob
import pysam
import sys

sample_l = []
sex_d = {}

# 1. 수정된 부분: 하위 디렉토리(sample_name) 안에 있는 dedup.bam 파일을 찾도록 패턴 변경
bams = glob.glob("*/*.dedup.bam")

print("\n************************************")
print("Sex determination script opening.")
print("************************************\n")

for i in bams:
    # i의 예시: "sample_name/sample_name.dedup.bam"
    # 2. 수정된 부분: 경로에서 파일명만 먼저 추출한 뒤 분리하는 것이 더 안전합니다.
    sample = os.path.basename(i).split(".")[0]
    
    sample_l.append(sample)
    bamfile = pysam.AlignmentFile(i, "rb", check_sq=False)
    sry_count = bamfile.count(contig='chrY', start=2786989, stop=2787603, until_eof=False, read_callback='all')
    
    # 어떤 샘플의 count인지 알기 쉽게 출력 포맷을 살짝 변경했습니다.
    print(f"[{sample}] SRY count: {sry_count}")

    if sry_count >= 50:
        sex = "M"
    elif sry_count <= 10:
        sex = "F"
    else:
        sex = "?"
    
    sex_d[sample] = sex

if os.path.isfile("samples.txt"):
    print("Sex determination was already done. No changes.")
else: 
    with open('samples.txt', 'w') as samp:
        samp.write("sample\tsex\n")
        
        for s in sample_l:
            samp.write(s + "\t" + sex_d[s] + "\n")
        
        print("Sex determination done for {} samples.".format(str(len(sample_l))))
        print("samples.txt generated.")

print("\nSex determination script job done!\n")
