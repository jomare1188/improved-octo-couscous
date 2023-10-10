rm summary.txt
grep "#Reads" *.refstats|sed -r 's/:#Reads\t/\tInputReads\t/' | sed 's/.refstats//' > summary.txt
grep "#Mapped" *.refstats|sed -r 's/:#Mapped\t/\tMapped\t/' | sed 's/.refstats//' >> summary.txt
grep -v "#" *.refstats | grep -v 'adapters      '| cut -f 1,6  | sed 's/.refstats:/\t/' >> summary.txt
