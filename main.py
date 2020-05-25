# -*- coding: utf-8 -*-
import Smith_Waterman as sw
import random
import os

print("start")
refFile = open("ref.txt", "r")
readFile = open("read.txt", "r")
refSeq = refFile.read()
readSeq = readFile.read()
print("files opened")

refSeq = refSeq.replace('\n','')
readSeq = readSeq.replace('\n','')
print("car return removed\n")

swInstance=sw.Smith_Waterman(refSeq, print_details=True, pos_etudiees={}, get_align_seq=False)
swInstance.align(readSeq)

print(swInstance)

refFile.close()
readFile.close()