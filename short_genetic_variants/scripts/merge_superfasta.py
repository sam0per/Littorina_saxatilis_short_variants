import glob
import fileinput
import sys

end = 'N' * 900

docs = glob.glob('genome/NONtarget.fasta')
for name in docs:
    FI = open(name, 'r')
    FO = open(name.replace('fasta', 'fa'), 'w')
    for num, line in enumerate(FI):
    	if not line.startswith('>'):
        	
    		FO.write(line.rstrip('\n') + end)

    FI.close()
    FO.close()

	