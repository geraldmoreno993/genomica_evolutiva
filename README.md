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



# Codigo 6: Visualización en Artemis parte II
```
# 6.1: para obervar los headers de cada contig en todos los archivos *.fa
grep ">" *.fa

# 6.2: para observar todo el contenido de todos los archivos *.fa
cat *.fa

# 6.3: para obervar las 10 primeras lineas de cada archivo *.fa
cat OQ603638.fa | head -n 10
cat OQ603651.fa | head -n 10
cat SRR30716253.fa | head -n 10
cat SRR30716253.fa | head -n 10

# 6.4: para observar los resultados de la anotación en ARTEMIS, debe contar con el archivo *fa original y el  archivo *.gff. Cargue primero el genoma y luego el archivo de anotación
conda activate art
art
file >> open file manager >> cargar el genoma en extension *fa
file >> read and entry >> cargar el archivo *gff
explorar

# 6.5: identifique las regiones inferidas por el programa (CDS), identifique si se identificó la identidad de esas regiones o si algunas aparecen como "hypothetical"
```

# codigo 7 : Ensamblaje Nanopore (programas)


```
# 7.1 : instalar los programas
# 7.1.1 : NanoPlot : calidad de secuencias Nanopore

conda install -c conda-forge -c bioconda nanoplot

or

pip install NanoPlot
pip install NanoPlot --upgrade

# 7.1.2 : Nanofilt : Filtrado por calidad de lecturas Nanopore

conda install -c bioconda nanofilt

or 

pip install nanofilt
pip install nanofilt --upgrade

# 7.1.3 : Flye: de-novo assembly

conda install -c bioconda flye

or 

git clone https://github.com/fenderglass/Flye
cd Flye
python setup.py install

# 7.1.4 : Minimap2 : polishing (parte 1)

conda install -c bioconda minimap2

or

git clone https://github.com/lh3/minimap2
cd minimap2 && make

# 7.1.5 : Racon : polishing (parte 2)

conda install -c bioconda racon

or 

git clone --recursive https://github.com/lbcb-sci/racon.git racon
cd racon
mkdir build
cd build
cmake -DCMAKE_BUILD_TYPE=Release ..
make
cd build/bin/ 
export PATH=$PATH:$HOME/bin
cp racon $HOME/bin
chmod +x $HOME/bin/racon

# 7.1.6 : Requerimientos de MEDAKA (Pyabpoa, bcftools, samtools (v1.11), minimap2)

pip install pyabpoa
sudo apt install bcftools
conda install -c bioconda samtools==1.11

# 7.1.7 : MEDAKA, secuencias consenso (si MEDAKA no funciona correctamente, instala los programas requeridos)

conda install -c conda-forge –c bioconda medaka

or

pip install medaka
```

# codigo 8 : ensamblaje Nanopore (pipeline)
material de apoyo > https://denbi-nanopore-training-course.readthedocs.io/en/stable/index.html
```
# 8.1: descargar la informacion (códigos SRR17110067 y SRR17110070)
mkdir sra_files ;
prefetch --max-size 50G --option-file accessions.txt ;
mv */*.sra . ;
fasterq-dump --split-files *.sra 
gzip *.fastq ;
mkdir sra_files ;
mv *.sra sra_files/ ;

# 8.2: inspeccionar las longitudes de los reads ##
zcat SRR17110067.fastq.gz | grep -n "length" | cut -f2 -d'=' | sort -r -n | uniq | head -n 20
zcat SRR17110070.fastq.gz | grep -n "length" | cut -f2 -d'=' | sort -r -n | uniq | head -n 20

# 8.3: NanoPlot
NanoPlot -t 2 -o SRR17110067_QC --fastq SRR17110067.fastq.gz
NanoPlot -t 2 -o SRR17110070_QC --fastq SRR17110070.fastq.gz

# 8.4: NanoFilt
gunzip -c SRR17110067.fastq.gz | NanoFilt --logfile nanofilt.log -l 500 -q 10 | gzip > SRR17110067.trim.fastq.gz ;
gunzip -c SRR17110070.fastq.gz | NanoFilt --logfile nanofilt.log -l 500 -q 10 | gzip > SRR17110070.trim.fastq.gz ;
ls -lh ;

# 8.5: Flye
flye -o SRR17110067.genoma --nano-raw SRR17110067.trim.fastq.gz --threads 4 ;
flye -o SRR17110070.genoma --nano-raw SRR17110070.trim.fastq.gz --threads 4 ;
ls -lh ;

# 8.6 : Minimap2 + Racon (Polishing)
minimap2 -x ava-ont -t 4 SRR17110067.genoma/assembly.fasta SRR17110067.trim.fastq.gz > overlaps1.paf ;
racon -t 4 SRR17110067.trim.fastq.gz overlaps1.paf SRR17110067.genoma/assembly.fasta > SRR17110067.racon1.fasta ;

minimap2 -x ava-ont -t 4 SRR17110070.genoma/assembly.fasta SRR17110070.trim.fastq.gz > overlaps2.paf ;
racon -t 4 SRR17110070.trim.fastq.gz overlaps2.paf SRR17110070.genoma/assembly.fasta > SRR17110070.racon1.fasta ;

minimap2 -x ava-ont -t 4 SRR17110067.racon1.fasta SRR17110067.trim.fastq.gz > overlaps3.paf ;
racon -t 4 SRR17110067.trim.fastq.gz overlaps3.paf SRR17110067.racon1.fasta > SRR17110067.racon2.fasta ;

minimap2 -x ava-ont -t 4 SRR17110070.racon1.fasta SRR17110070.trim.fastq.gz > overlaps4.paf ;
racon -t 4 SRR17110070.trim.fastq.gz overlaps4.paf SRR17110070.racon1.fasta > SRR17110070.racon2.fasta ;

# 8.7 : Medaka (consensus)
medaka_consensus -i SRR17110070.trim.fastq.gz -d SRR17110070.racon2.fasta -o medaka_SRR17110070 -t 4 ;
medaka_consensus -i SRR17110067.trim.fastq.gz -d SRR17110067.racon2.fasta -o medaka_SRR17110067 -t 4 ;

# 8.8 : QUAST
quast.py -o quast_results -m 0 consensus.fasta

# 8.9 : Bandage
```


# codigo 9 : BLAST (pipeline)
material de apoyo > (http://www.mgc.ac.cn/VFs/download.htm)
Lo que se hará es construir una base de datos local con todos los factores de virulencia (nucleotidos y aminoacidos) y luego usar como query los genomas de mycoplasma para 
blastear factores de virulencia
```
pwd

nano Chlamydia_accessions.txt

wget https://ftp.ncbi.nlm.nih.gov/pub/datasets/command-line/v2/linux-amd64/datasets
wget https://ftp.ncbi.nlm.nih.gov/pub/datasets/command-line/v2/linux-amd64/dataformat
./datasets lit Chlamydia_accessions.txt
# 9.1 : instalacion a traves de CONDA
conda create --name blast
conda install bioconda::blast



# Descarga en formato GenBank y guarda en un archivo zip
./datasets download genome accession --inputfile Chlamydia_accessions.txt --filename Chlamydia_fasta.zip --include genome,seq-report
unzip Chlamydia_fasta.zip -d Chlamydia_data/

# Alternativamente, descarga en formato FASTA
./dataformat tsv genome --package Chlamydia_fasta.zip --fields organism-name,assminfo-level,accession > Chlamydia_metadata.tsv


mv ./*/*.fna .

#Descarga de la totalidad de la base de datos vfdb
wget http://www.mgc.ac.cn/VFs/Down/VFDB_setB_nt.fas.gz
wget http://www.mgc.ac.cn/VFs/Down/VFDB_setB_pro.fas.gz


conda activate blast

conda install bioconda::blast

gunzip *
#RUN BLAST

# 9.7 : run BLAST+
makeblastdb -in VFDB_setB_nt.fas -dbtype nucl ;
blastn -db VFDB_setB_nt.fas -query GCA_001183825.1_ASM118382v1_genomic.fna -perc_identity 90 -outfmt 6 -num_threads 4 > blast.csv ; ##El mejor formato es el 10, (csv), el 6 es el (tab)
head blast.csv ;
cat blast.csv ;

# 9.8 : headers (insertar un header)
sed '1i query.acc.ver subject.acc.ver perc.identity alignment.length mismatches gap.opens q.start q.end s.start s.end evalue bit.score' blast.csv | tr " " "\t" > blast.2.csv

# 9.9 : revisar resultados
head blast.2.csv
cat blast.2.csv
wc -l blast.2.csv
```
