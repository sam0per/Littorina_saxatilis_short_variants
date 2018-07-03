import glob
import fileinput
import sys


outs = glob.glob('genome/NONtarget.fa')
for num, line in enumerate(fileinput.input(outs, inplace=True)):
	sys.stdout.write('>superscaff' + str(num) + '\n{l}'.format(l=line))



