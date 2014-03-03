import sys
file_name=sys.argv[1]

with open(file_name+'.new','w') as output_file:
    with open(file_name,'r') as read_file:
        for line in read_file:
            info=line.split('\t')
            if info[4].strip()=='Pathogenic':
                info[4]='1\n'
            else:
                info[4]='0\n'
            new_line=''
            for i in range(0,4):
                new_line=new_line+info[i]+'\t'
            new_line=new_line+info[4]
            output_file.write(new_line)
