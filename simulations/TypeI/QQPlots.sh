source $1
grep $gene $null | cut -f2 > /tmp/1.txt
grep $gene $recomb | cut -f2 > /tmp/2.txt
grep $gene $ldmissing | cut -f2 > /tmp/3.txt
echo """ 
d1 = sort(scan('/tmp/1.txt'))
d2 = sort(scan('/tmp/2.txt'))
d3 = sort(scan('/tmp/3.txt'))
d20 = sort(scan('/tmp/1.txt')[1:length(d2)])
d30 = sort(scan('/tmp/1.txt')[1:length(d3)])
print(length(which(d1>3.6)))
pdf('$gene.pdf')
hist(d1)
plot(d20,d2)
abline(0,1)
plot(d30,d3)
abline(0,1)
dev.off()
""" | R --slave  
if [ "$(uname)" == "Darwin" ]; then
    prog='open'
else
    prog='google-chrome'
fi
$prog $gene.pdf 
