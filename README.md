# genomica_evolutiva
Clases de genómica evolutiva de la maestría en Bioinformática y Ciencias Ómicas UNMSM 2024


#Código 1: Descargar sratoolkit

```
mkdir genomas;
grep ">" genomas.fasta ;
ls ;

wget https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/3.1.1/sratoolkit.3.1.1-ubuntu64.tar.gz -O stk.tar.gz
chmod +x stk.tar.gz #activar permisos
tar -vxzf stk.tar.gz  #comprimir

export PATH=$PATH:$PWD/sratoolkit.3.0.10-ubuntu64/bin

```


#Código 2: Descargar archivo .sra y dividirlo en dos archivos .fastq

```
prefetch -h 
prefetch --max-size 50G --option-file sra_accessions_1.txt  #;
rm -r ERR12389866/ ERR12543675/  ##Borrar las carpetas pero con cuidado
fasterq-dump --split-files *.sra #separar en forward y reverse;
gzip *fastq ;
fastqc *
```
