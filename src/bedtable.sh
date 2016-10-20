
bedtable.head(){
usage="$FUNCNAME <bedtable> [option]
[options]: 
	-v : output body
";
if [ $# -lt 1 ];then echo "$usage"; return; fi
	cat $1 | awk -v OFS="\t" -v O=$2 '{
		if( $1 ~ /^#/){ # do nothing 
		}else{
			b=0;
			if( NF > 5 && $2~/^[0-9]+$/){ b=1;}

			if(O == "-v"){
				if( b==1) print $0;
			}else{
				if(b==0) print $0;
				else exit 0;
			}
		}
	}'
}
bedtable.head.test(){
echo \
"# hi
chrom	start	end	name	score	strand	score1	score2
chr1	1	2	n1	0	+	1	1
chr1	1	2	n1	0	+	1	1" > tmp.a
echo "head:"
bedtable.head tmp.a 
echo "body:"
bedtable.head tmp.a -v
rm tmp.*
}


bedtable.intersect(){
usage="$FUNCNAME <a.bedtable> <b.bedtable> <output> [options]
[options]: same as the options in intersectBed 
";
opts="${@:4}"
if [ $# -lt 3 ];then echo "$usage"; return; fi
	bedtable.head $1 -v > $3.abed
	bedtable.head $2 -v > $3.bbed
	intersectBed -a $3.abed -b $3.bbed $opts  > $3
	rm $3.*	
}

bedtable.test(){
echo \
"# hi
chrom	start	end	name	score	strand	score1	score2
chr1	1	2	n1	0	+	1	1
chr1	1	2	n1	0	+	1	1
chr1	1	2	n1	0	+	1	1" > tmp.a
bedtable.head tmp.a -v > tmp.b
bedtable.intersect tmp.a tmp.a tmp.o -wa -wb -s 
head tmp.o
}

