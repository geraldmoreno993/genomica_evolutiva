# GENOMICA EVOLUTIVA
Clases de genómica evolutiva de la maestría en Bioinformática y Ciencias Ómicas UNMSM 2024


# Código 1: Descargar sratoolkit

```
mkdir genomas;
grep ">" genomas.fasta ;
ls ;

wget https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/3.1.1/sratoolkit.3.1.1-ubuntu64.tar.gz -O stk.tar.gz
chmod +x stk.tar.gz #activar permisos
tar -vxzf stk.tar.gz  #descomprimir

export PATH=$PATH:$PWD/sratoolkit.3.0.10-ubuntu64/bin

```


# Código 2: Descargar archivo .sra y dividirlo en dos archivos .fastq

```
prefetch -h 
prefetch --max-size 50G --option-file sra_accessions_1.txt  #;
rm -r ERR12389866/ ERR12543675/  ##Borrar las carpetas pero con cuidado
fasterq-dump --split-files *.sra #separar en forward y reverse;
gzip *fastq ;
fastqc *
```


# Código 3: Ensamblaje por mapeo
## Clase de 23 de set 2024

```


##Tratar de que todos los archivos descargados Fastq, referencia, sam y los indexados esten dentro de una misma carpeta
#0# descargar de NCBI
prefetch --max-size 50G --option-file accessions_mpox.txt ;
mv */*.sra . ;
fasterq-dump --split-files *.sra ;
gzip *fastq ;
fastqc * ;

#1# indexar el genoma de referencia#
bwa index reference.fasta ;

##Ejecutar todo el paso 2 y 3 en una sola corrida porqeu es un bucle para obtener archivo .sam y bam ########################
#2# preparar las instrucciones generales#
for r1 in *fastq.gz
do
prefix=$(basename $r1 _1.fastq.gz)
r2=${prefix}_2.fastq.gz

#3# instrucciones para generar el archivo .bam#
bwa mem -t 4 reference.fasta $r1 $r2 > ${prefix}_uno.sam ;
samtools view -@ 4 -bS -T reference.fasta ${prefix}_uno.sam > ${prefix}_unoa.bam ;
samtools sort -@ 4 -n ${prefix}_unoa.bam -o ${prefix}_dosa.bam ;
samtools fixmate -@ 4 -m ${prefix}_dosa.bam ${prefix}_tresa.bam ;
samtools sort -@ 4 ${prefix}_tresa.bam -o ${prefix}_cuatroa.bam ;
samtools markdup -@ 4 ${prefix}_cuatroa.bam ${prefix}.bam ;
samtools index -@ 4 ${prefix}.bam ;
rm ${prefix}_uno.sam ${prefix}_unoa.bam ${prefix}_dosa.bam ${prefix}_tresa.bam ${prefix}_cuatroa.bam ;
done ;
ls ;
########################################

#4# extraer genoams consenso #
for r1 in *bam
do
prefix=$(basename $r1 .bam)
#2#estimate Ns#
samtools mpileup -aa -A -d 0 -Q 0 $r1 | ivar consensus -p ${prefix}.fasta -q 25 -t 0.6 -m 10 ;
done ; 
ls ;

#5# annotation #
mkdir -p annotation ;
mkdir -p ffn ;
for r1 in *fa
do
prefix=$(basename $r1 .fa)
prokka --cpus 4 $r1 -o ${prefix} --prefix ${prefix} --kingdom Viruses ; 
mv ${prefix}/*.gff annotation/${prefix}.gff
done ;
conda deactivate ;
cp */*.ffn ffn/ ; 
ls ; 



```

# Código 4: Reporte IGV y anotaciòn con Prokka
## Clase de 30 de set 2024


```
#######
##IGV##
#######

conda install bioconda::igv

igv
#Software con interfaz gráfica que permite colocar como input el genoma de referencia en formato fasta y el bam.bai de la secuencia que vas a alinear o emsamblar.

############
###PROKKA###
############

conda deactivate
conda activate prokka
cd /home/gerald/Documentos/maestria/2do_ciclo/genomica_evolutiva/clase4
ls

#Prokka generalmente acepta archivos.fasta, cambiando de extension a .fasta
for file in *.fa; do
    mv -- "$file" "${file%.fa}.fasta"
done


#Bucle para anotar multiples archivos fasta
for a1 in *.fasta
do
prefix=$(basename $a1 .fasta)
prokka --cpus 4 $a1 -o ${prefix} --prefix ${prefix} --kingdom Viruses ;
done ;

```


# Código 5: Visualización de genoma anotado con Artemis
## Clase de 3 de octubre 2024



```
# instalacion de artemis #
conda create -n art
conda activate art
conda install bioconda::artemis

#Invocar el programa con:

art

#Dentro del programa con interfaz gráfica cargar en:
#file > open file manager > cargar tu archivo de anotaicion .gff
#en la ventana inferior aparecen los CDS haz click derecho y coloca show gene names 



```

