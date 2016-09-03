
bed.enc(){
## encode bed6 in the name field 
	cat $1 | awk -v OFS="\t" -v j=${2:-"@"} '{ n=$1;for(i=2;i<= 6;i++){ n=n j $(i);} $4=n; }1' 
}


bed.3p(){
 	awk -v OFS="\t" '{ if($6=="-"){$3=$2+1;} $2=$3-1; print $0; }' $1
}

bed.5p(){
 	awk -v OFS="\t" '{if($6=="-"){$2=$3-1;}$3=$2+1; print $0; }' $1
}

bed.flank(){
usage="
FUNCT : extract flanking regions
USAGE : $FUNCNAME <bed> <left> <right> <strand_opt>
	<strand_opt> :  -s: strand specific
" 
if [ $# -ne 4 ];then echo "$usage"; return; fi
	awk -v OFS="\t" -v L=$2 -v R=$3 -v S=$4 '{ 
		if(S == "-s"  && $6 == "-"){
			$2=$2-R; $3=$3+L;
		}else{ $2=$2-L; $3=$3+R; } 
		if( $2 < 0){ $2=0;} if( $3 < 0){ $3=1;}
		print $0;
	}' $1;
}
