#!/usr/bin/env bash

Illumina_adapters=/home/amanda/miniconda3/share/trimmomatic-0.39-2/adapters/TruSeq3-PE-2.fa
run_interpro_sequence_dir=/home/amanda/PROYECTS/Mauri_Imparalogos/Interpro_Analisys/sequences_all
run_interpro_imparalogs_main_result=/home/amanda/PROYECTS/Mauri_Imparalogos/Interpro_Analisys/inparalogos_info
interprot_path=/home/amanda/PROYECTS/interproscan-5.56-89.0
interpro_ID_info=interpro.xml

schisto_genome_fasta=schistosoma_mansoni.PRJEA36577.WBPS18.genomic.fa
schisto_gff_file=schistosoma_mansoni.PRJEA36577.WBPS18.annotations.gff3
schisto_trans=schistosoma_mansoni.PRJEA36577.WBPS18.mRNA_transcripts.fa
eggnog_annotation=Eggnog_annotation.in

Schisto_TRIM=FALSE
Schisto_STAR=FALSE
Schisto_COUNT=FALSE
Schisto_TABLE_MAKER=FALSE

Schisto_Protasio_TRIM=FALSE
Schisto_Protasio_STAR=FALSE
Schisto_Protasio_COUNT=FALSE
Schisto_Protasio_TABLE_MAKER=FALSE

Schisto_Sanger_TRIM=FALSE
Schisto_Sanger_STAR=FALSE
Schisto_Sanger_COUNT=FALSE
Schisto_Sanger_TABLE_MAKER=FALSE

Smansoni_Genome_Location=FALSE
Smansoni_Destance_groups=FALSE


Dif_Exp_Table=FALSE ##### No utilizado
INTERPRO_PREPARE=FALSE
INTERPRO_RUN=FALSE
INTERPRO_DOMAINS_IMP=FALSE
INTERPRO_GO_TERMS=FALSE
INTERPRO_DOMAINS_PRE_HEATMAP=FALSE

SUP_Tables=FALSE
SUP_ANOTATION=FALSE
SUP_ANOTATION_CONSULTA=TRUE

confirmacion_recorrida=FALSE
Final_Summary=FALSE

threads=30
run_trim () {
  rm -r $species_name"/Clean_Reads"
  mkdir $species_name"/Clean_Reads"

  read_files_illumina=$(ls $species_name"/raw_reads/"*"_1.fastq.gz" | sed "s/.*\///" | sed "s/_1.fastq.gz//")
  echo $species_name" Trimming Reads"
  for read in $read_files_illumina
  do
    echo "Procesando: "$read
    trimmomatic PE $species_name"/raw_reads/"$read"_1.fastq.gz" $species_name"/raw_reads/"$read"_2.fastq.gz" -baseout $species_name"/Clean_Reads/"$read'_trimmed_read.gz' -threads $threads ILLUMINACLIP:$Illumina_adapters":2:30:10" LEADING:3 TRAILING:3 MINLEN:36
  done
  seqkit stats -a -T $species_name"/Clean_Reads/"*".gz" >> $species_name"/Clean_Reads/Basic_metrics.tab"
}

con_multimaping () {
  # Haciendo referencias y mapeos
  rm -r $species_name"/Multimaping"
  mkdir $species_name"/Multimaping"
  rm -r $species_name"/Multimaping/Counts"
  mkdir $species_name"/Multimaping/Counts"
  mkdir $species_name"/Multimaping/STAR_aln"
  mkdir $species_name"/Multimaping/STAR_aln/Referencias"

  read_files_illumina=$(ls $species_name"/Clean_Reads/"*"_trimmed_read_1P.gz" | sed "s/.*\///" | sed "s/_trimmed_read_1P.gz//")
  STAR --runMode genomeGenerate --runThreadN $threads --genomeSAindexNbases 12 --genomeDir $species_name"/Multimaping/STAR_aln/Referencias" --genomeFastaFiles $genome_fasta --sjdbGTFfile $gff_file
  for read in $read_files_illumina
  do
    echo "Descomprimiendo"

    echo "Preparando papaya para STAR: "$read
    gzip -d -k $species_name"/Clean_Reads/"$read'_trimmed_read_1P.gz'

    echo "Preparando papaya para STAR: "$read
    gzip -d -k $species_name"/Clean_Reads/"$read'_trimmed_read_2P.gz'

    echo "Mapeando"
    STAR --runMode alignReads --runThreadN $threads --genomeSAindexNbases 12 --outFileNamePrefix $species_name"/Multimaping/STAR_aln/"$read"_aln" --outSAMtype BAM SortedByCoordinate --genomeDir  $species_name"/Multimaping/STAR_aln/Referencias"  --sjdbGTFfile $gff_file --limitBAMsortRAM 6213894447 --readFilesIn $species_name"/Clean_Reads/"$read'_trimmed_read_1P' $species_name"/Clean_Reads/"$read'_trimmed_read_2P'

    rm $species_name"/Clean_Reads/"$read'_trimmed_read_1P'
    rm $species_name"/Clean_Reads/"$read'_trimmed_read_2P'

    samtools index $species_name"/Multimaping/STAR_aln/"$read"_alnAligned.sortedByCoord.out.bam"
    htseq-count -r pos -s no --idattr=gene_id --format=bam --type=exon -q $species_name"/Multimaping/STAR_aln/"$read"_alnAligned.sortedByCoord.out.bam" $species_name"/temp_"$species_name".gtf" >> $species_name"/Multimaping/Counts/"$read".count"
  done

  extract_count=$(awk -F "\t" '{print $1}' $species_name"/Multimaping/Counts/"*".count"| sort -u | grep -v __)
  echo "Term_Count "$read_files_illumina | tr " " "\t" >> $species_name"/Multimaping/Counts/"$species_name"_all_counts.tab"
  for term in $extract_count
  do
    echo "Making: "$term
    registro=$(echo $term | sed "s/gene://")
    for read in $read_files_illumina
    do
      extract=$(grep -F -w $term $species_name"/Multimaping/Counts/"$read".count" | awk -F "\t" '{print $2}')
      registro=$(echo $registro" "$extract)
    done
    echo $registro | tr " " "\t" >> $species_name"/Multimaping/Counts/"$species_name"_all_counts.tab"
  done
}

star_aln () {
  # Haciendo referencias y mapeos
  rm -r $species_name"/STAR_aln"
  mkdir $species_name"/STAR_aln"
  mkdir $species_name"/STAR_aln/Referencias"

  read_files_illumina=$(ls $species_name"/Clean_Reads/"*"_trimmed_read_1P.gz" | sed "s/.*\///" | sed "s/_trimmed_read_1P.gz//")
  STAR --runMode genomeGenerate --runThreadN $threads --genomeSAindexNbases 12 --genomeDir $species_name"/STAR_aln/Referencias" --genomeFastaFiles $genome_fasta --sjdbGTFfile $gff_file
  for read in $read_files_illumina
  do
    echo "Descomprimiendo"

    echo "Preparando papaya para STAR: "$read
    gzip -d -k $species_name"/Clean_Reads/"$read'_trimmed_read_1P.gz'

    echo "Preparando papaya para STAR: "$read
    gzip -d -k $species_name"/Clean_Reads/"$read'_trimmed_read_2P.gz'

    echo "Mapeando"
    # STAR --runMode alignReads --runThreadN $threads --genomeSAindexNbases 12 --outFileNamePrefix $species_name"/STAR_aln/"$read"_aln" --outSAMtype BAM SortedByCoordinate --genomeDir  $species_name"/STAR_aln/Referencias" --outFilterMultimapNmax 1 --sjdbGTFfile $gff_file --readFilesIn $species_name"/Clean_Reads/"$read'_trimmed_read_1P' $species_name"/Clean_Reads/"$read'_trimmed_read_2P'
    STAR --runMode alignReads --runThreadN $threads --genomeSAindexNbases 12 --outFileNamePrefix $species_name"/STAR_aln/"$read"_aln" --outSAMtype BAM SortedByCoordinate --genomeDir  $species_name"/STAR_aln/Referencias" --outFilterMultimapNmax 1 --sjdbGTFfile $gff_file --limitBAMsortRAM 6213894447 --readFilesIn $species_name"/Clean_Reads/"$read'_trimmed_read_1P' $species_name"/Clean_Reads/"$read'_trimmed_read_2P'

    rm $species_name"/Clean_Reads/"$read'_trimmed_read_1P'
    rm $species_name"/Clean_Reads/"$read'_trimmed_read_2P'

    samtools index $species_name"/STAR_aln/"$read"_alnAligned.sortedByCoord.out.bam"
  done
}

gff_shenanigas () {
  rm $species_name"/"$species_name"_tran_len.tab" $species_name"/temp_"$species_name".gtf" $species_name"/work_"$species_name".gtf"
  echo "Seqkit!!"
  seqkit fx2tab -n -i -l $mrna >> $species_name"/"$species_name"_tran_len.tab"
  echo "GFFread!!!"
  gffread $schisto_gff_file -F -C -T -o $species_name"/temp_"$species_name".gtf"
  echo "Done GFFread"

  genes_id=$(awk -F "\t" '{if ($3 == "transcript") print $9}' $species_name"/temp_"$species_name".gtf" | tr ";" "\n" | grep gene_id | awk -F ":" '{print $2}' | tr -d '"' | sort -u)
  for gene in $genes_id
  do
    check_num=$(grep "gene:"$gene\"\; $species_name"/temp_"$species_name".gtf" | awk -F "\t" '{if ($3 == "transcript") print $9}' | grep -c .)
    echo $gene" "$check_num
    if [ $check_num -eq 1 ]
    then
      grep -w "gene:"$gene\"\; $species_name"/temp_"$species_name".gtf" >> $species_name"/work_"$species_name".gtf"
      grep $gene $species_name"/"$species_name"_tran_len.tab" >> $species_name"/"$species_name"_tran_len_best.tab"
    else
      best_trans=$(grep $gene $species_name"/"$species_name"_tran_len.tab" | sort -n -k 2 | tail -n 1 | awk -F "\t" '{print $1}')
      echo "this happens!!" $best_trans
      grep -w "transcript:"$best_trans\"\; $species_name"/temp_"$species_name".gtf" >> $species_name"/work_"$species_name".gtf"
      grep -w $best_trans $species_name"/"$species_name"_tran_len.tab" >> $species_name"/"$species_name"_tran_len_best.tab"
    fi
  done
}

counting () {
  rm -r $species_name"/Counts"
  mkdir $species_name"/Counts"

  read_files_illumina=$(ls $species_name"/Clean_Reads/"*"_trimmed_read_1P.gz" | sed "s/.*\///" | sed "s/_trimmed_read_1P.gz//")
  for read in $read_files_illumina
  do
    echo "Contando Mapeo: "$species_name"/STAR_aln/"$read"_alnAligned.sortedByCoord.out.bam"
    echo "GFF: "$species_name"/temp_"$species_name".gtf"
    echo ""
    htseq-count -r pos -s no --idattr=gene_id --format=bam --type=exon -q $species_name"/STAR_aln/"$read"_alnAligned.sortedByCoord.out.bam" $species_name"/work_"$species_name".gtf" >> $species_name"/Counts/"$read".count"
  done
}

table_make () {
  rm -r $species_name"/Work_Tables"
  mkdir $species_name"/Work_Tables"

  read_files_illumina=$(ls $species_name"/Clean_Reads/"*"_trimmed_read_1P.gz" | sed "s/.*\///" | sed "s/_trimmed_read_1P.gz//")
  extract_count=$(awk -F "\t" '{print $1}' $species_name"/Counts/"*".count" | sort -u | grep -v __)

  echo "Term_Count "$read_files_illumina | tr " " "\t" >> $species_name"/Work_Tables/"$species_name"_all_counts.tab"
  for term in $extract_count
  do
    echo "Making: "$term
    registro=$(echo $term | sed "s/gene://")
    for sample in $read_files_illumina
    do
      extract=$(grep -F -w $term $species_name"/Counts/"$sample".count" | awk -F "\t" '{print $2}')
      registro=$(echo $registro" "$extract)
    done
    echo $registro | tr " " "\t" >> $species_name"/Work_Tables/"$species_name"_all_counts.tab"
  done
}

genome_context () {
  rm -r $species_name"/Genome_Locations"
  mkdir $species_name"/Genome_Locations"
  mkdir $species_name"/Genome_Locations/Temporal"

  count=1
  recorrer=$(grep -c . $imparalog_file)

  while [ $count -le $recorrer ]
  do
    frak_you_mauri=
    imp_ID=$(sed -n $count"p" $imparalog_file | awk -F "\t" '{print $1}')
    genes=$(sed -n $count"p" $imparalog_file | awk -F "\t" '{print $2}' | tr "," "\n" | tr -d " " | awk -F "." '{print $1}')
    total_genes=$(echo $genes | tr " " "\n" | grep -c .)
    total_safe_genes=$(echo $genes | tr " " "\n" | grep -c .)

    echo $imp_ID": "$count"/"$recorrer

    echo "Imparalog ID: "$imp_ID" ("$imp_ID"WORMBASE_UPDATE)" >> $species_name"/Genome_Locations/Genomic_Context.txt"
    echo "Genes:"$genes | sed "s/ /; /g" | sed "s/:/: /" >> $species_name"/Genome_Locations/Genomic_Context.txt"
    echo "Status:"$imp_ID"XXXCHANGE_MEXXXX" | sed "s/ /; /g" | sed "s/:/: /" >> $species_name"/Genome_Locations/Genomic_Context.txt"
    echo "Not part of new update: "$imp_ID"YOU_FRAKING_TOASTER" >> $species_name"/Genome_Locations/Genomic_Context.txt"

    echo "" >> $species_name"/Genome_Locations/Genomic_Context.txt"

    for gen in $genes
    do
      echo "ID=gene:"$gen";" >> $species_name"/Genome_Locations/Temporal/Exclude.temp"
      grep "ID=gene:"$gen";" $gff_file | awk -F "\t" -v gen=$gen '{if ($3=="gene") print gen"\t"$1"\t"$4"\t"$5"\t"$7"\tImparalogo"}' | grep -v "#" >> $species_name"/Genome_Locations/Temporal/Coordinates.temp"

      check_update=$(grep -w -F -c $gen $new_annotation)
      if [ $check_update -eq 0 ]
      then
        total_safe_genes=$(($total_safe_genes - 1))
        frak_you_mauri=$(echo $frak_you_mauri";"$gen | sed "s/^;//")
        echo $imp_ID" "$gen >> $species_name"/Genome_Locations/Update_warnings.txt"
      fi
    done

    imp_status=
    vibe_check="Good"
    target_chr=$(awk -F "\t" '{print $2}' $species_name"/Genome_Locations/Temporal/Coordinates.temp" | sort -u)
    for chr in $target_chr
    do
      grep -w -F $chr $species_name"/Genome_Locations/Temporal/Coordinates.temp" | grep -v "#" >> $species_name"/Genome_Locations/Temporal/Coordinates2.temp"
      start_search=$(grep -w -F $chr $species_name"/Genome_Locations/Temporal/Coordinates.temp" | sort -n -k 3 | head -n 1 | awk -F "\t" '{print $3}')
      end_search=$(grep -w -F $chr $species_name"/Genome_Locations/Temporal/Coordinates.temp" | sort -n -k 3 | tail -n 1 | awk -F "\t" '{print $4}')

      gffread -O -R -r $chr":"$start_search".."$end_search $gff_file | awk -F "\t" '{if ($3!="CDS" && $3!="exon" && $3!="CDS" && $3!="mRNA") print }' | sed "s/;.*//" | sed "s/ID=gene://" | awk -F "\t" '{print $9" "$1" "$4" "$5" "$7" extra_extra"}' | grep -w -v -f $species_name"/Genome_Locations/Temporal/Exclude.temp" | grep -v "#" | tr " " "\t" >> $species_name"/Genome_Locations/Temporal/Coordinates2.temp"

      grep -w -F $chr $species_name"/Genome_Locations/Temporal/Coordinates2.temp" | sort -n -k 3 >> $species_name"/Genome_Locations/Temporal/Coordinates3.temp"
      element_check=$(grep -w -c extra_extra $species_name"/Genome_Locations/Temporal/Coordinates3.temp")

      if [ $element_check -le 10 ]
      then
        int_count=1
        int_stop=$(grep -c . $species_name"/Genome_Locations/Temporal/Coordinates3.temp")
        registro_genes=
        registro="=="
        while [ $int_count -le $int_stop ]
        do
          int_gen=$(sed -n $int_count"p" $species_name"/Genome_Locations/Temporal/Coordinates3.temp" | awk -F "\t" '{print $1}')
          type=$(sed -n $int_count"p" $species_name"/Genome_Locations/Temporal/Coordinates3.temp" | awk -F "\t" '{print $6}')

          if [ $type == "extra_extra" ]
          then
            int_gen=$(echo $int_gen"*")
          fi
          registro_genes=$(echo $registro_genes";"$int_gen | sed "s/^;//")

          int_start=$(sed -n $int_count"p" $species_name"/Genome_Locations/Temporal/Coordinates3.temp" | awk -F "\t" '{print $3}')
          int_strand=$(sed -n $int_count"p" $species_name"/Genome_Locations/Temporal/Coordinates3.temp" | awk -F "\t" '{print $5}')

          next=$(($int_count + 1))
          if [ $next -le $int_stop ]
          then
            int_end=$(sed -n $next"p" $species_name"/Genome_Locations/Temporal/Coordinates3.temp" | awk -F "\t" '{print $4}')
            separation=$(($int_end - $int_start + 1))
            registro=$(echo $registro"=="$int_gen"["$int_strand"]==("$separation"bp)")
          else
            registro=$(echo $registro"=="$int_gen"["$int_strand"]====")
          fi
          int_count=$(($int_count + 1))
        done
        echo "------------------------------------------------" >> $species_name"/Genome_Locations/Genomic_Context.txt"
        echo $chr" ("$registro_genes")" >> $species_name"/Genome_Locations/Genomic_Context.txt"
        echo "" >> $species_name"/Genome_Locations/Genomic_Context.txt"
        echo $registro >> $species_name"/Genome_Locations/Genomic_Context.txt"
        echo "" >> $species_name"/Genome_Locations/Genomic_Context.txt"
      else
        echo "------------------------------------------------" >> $species_name"/Genome_Locations/Genomic_Context.txt"
        vibe_check="Bad"
        echo $chr" (too many hits)" >> $species_name"/Genome_Locations/Genomic_Context.txt"
        echo "" >> $species_name"/Genome_Locations/Genomic_Context.txt"
        echo "Not_Shown" >> $species_name"/Genome_Locations/Genomic_Context.txt"
        echo "" >> $species_name"/Genome_Locations/Genomic_Context.txt"
      fi
      echo "Gen Start End Strand" | tr " " "\t" >> $species_name"/Genome_Locations/Genomic_Context.txt"
      grep -w Imparalogo $species_name"/Genome_Locations/Temporal/Coordinates3.temp" | awk -F "\t" '{print $1"\t"$3"\t"$4"\t"$5}' >> $species_name"/Genome_Locations/Genomic_Context.txt"

      rm $species_name"/Genome_Locations/Temporal/Coordinates2.temp"
      rm $species_name"/Genome_Locations/Temporal/Coordinates3.temp"
    done
    total_chromosomes=$(echo $target_chr | tr " " "\n" | grep -c .)

    if [ $total_chromosomes -eq 1 ]
    then
      if [ $vibe_check == "Good" ]
      then
        imp_status=Candidate
      else
        imp_status=Complex
      fi
    else
      imp_status=Multiple_Chr
    fi
    replace="Status: "$imp_ID"XXXCHANGE_MEXXXX"
    new="Status: "$imp_status

    sed -i "s/$replace/$new/" $species_name"/Genome_Locations/Genomic_Context.txt"

    replace=$imp_ID"WORMBASE_UPDATE"
    new=$(echo $total_safe_genes" of "$total_genes)
    sed -i "s/$replace/$new/" $species_name"/Genome_Locations/Genomic_Context.txt"

    replace=$imp_ID"YOU_FRAKING_TOASTER"
    new=$frak_you_mauri

    sed -i "s/$replace/$new/" $species_name"/Genome_Locations/Genomic_Context.txt"

    echo "" >> $species_name"/Genome_Locations/Genomic_Context.txt"
    echo "#################################################################################" >> $species_name"/Genome_Locations/Genomic_Context.txt"
    echo "" >> $species_name"/Genome_Locations/Genomic_Context.txt"

    rm $species_name"/Genome_Locations/Temporal/"*".temp"
    count=$(($count + 1))
  done
}

correlations_by_proximity () {
  rm -r $species_name"/Correlations_By_proximity"
  mkdir $species_name"/Correlations_By_proximity"
  mkdir $species_name"/Correlations_By_proximity/Temporal"

  count=1
  recorrer=$(grep -c . $imparalog_file)

  while [ $count -le $recorrer ]
  do
    imp_ID=$(sed -n $count"p" $imparalog_file | awk -F "\t" '{print $1}')
    genes=$(sed -n $count"p" $imparalog_file | awk -F "\t" '{print $2}' | tr "," "\n" | tr -d " " | awk -F "." '{print $1}')

    registro=Done

    for gen1 in $genes
    do
      for gen2 in $genes
      do
        if [ ! "$gen1" == "$gen2" ]
        then
          echo $imp_ID": "$gen1" "$gen2
          current=$(echo $gen1" "$gen2 | tr " " "\n" | sort | tr "\n" "_" | sed "s/_$//")
          checking_done=$(echo $registro | grep -c $current)

          if [ $checking_done -eq 0 ]
          then
            registro=$(echo $registro" "$current)
            grep "ID=gene:"$gen1";" $gff_file | awk -F "\t" -v gen=$gen1 '{if ($3=="gene") print gen"\t"$1"\t"$4"\t"$5}' >> $species_name"/Correlations_By_proximity/Temporal/"$imp_ID"_"$gen1"_"$gen2"_gene_locations.txt"
            grep "ID=gene:"$gen2";" $gff_file | awk -F "\t" -v gen=$gen2 '{if ($3=="gene") print gen"\t"$1"\t"$4"\t"$5}' >> $species_name"/Correlations_By_proximity/Temporal/"$imp_ID"_"$gen1"_"$gen2"_gene_locations.txt"

            chr1=$(grep -w -F $gen1 $species_name"/Correlations_By_proximity/Temporal/"$imp_ID"_"$gen1"_"$gen2"_gene_locations.txt" | awk -F "\t" '{print $2}')
            chr2=$(grep -w -F $gen2 $species_name"/Correlations_By_proximity/Temporal/"$imp_ID"_"$gen1"_"$gen2"_gene_locations.txt" | awk -F "\t" '{print $2}')

            if [ "$chr1" == "$chr2" ]
            then
              start=$(sort -k 3 -n $species_name"/Correlations_By_proximity/Temporal/"$imp_ID"_"$gen1"_"$gen2"_gene_locations.txt" | awk -F "\t" '{print $4}' | head -n 1 )
              end=$(sort -k 3 -n $species_name"/Correlations_By_proximity/Temporal/"$imp_ID"_"$gen1"_"$gen2"_gene_locations.txt" | awk -F "\t" '{print $3}' | tail -n 1 )
              echo "Confirmar: "$imp_ID"_"$gen1"_"$gen2
              distancia=$(($end - $start))

              group=$(echo $end" "$start | awk -F " " '{if ($1-$2 <= 0) print "Other"; else if ($1-$2 <= 10000) print "Close"; else print "Far"}')
            else
              distancia=$chr1"/"$chr2
              group=Diferent_Chr
            fi
            # rm $species_name"/Correlations_By_proximity/Temporal/"$imp_ID"_"$gen1"_"$gen2"_gene_locations.txt"
            echo $imp_ID" "$gen1" "$gen2" "$distancia" "$group | tr " " "\t" >> $species_name"/Correlations_By_proximity/Gene_groups.txt"
          fi
        fi
      done
    done
    count=$(($count + 1))
  done
}

multimaping_compare (){
  count=1
  recorrer=$(grep -c . $imparalog_file)

  one_mapping_zero_count=0
  multimapping_mapping_zero_count=0
  all_missing=0
  no_multi_mappings=0
  weirdness=0

  echo "Resumen" > $species_name"/Multimaping/Counts/"$species_name"_comparisons.tab"
  echo "Total Sin datos (ambos): __REPLACE_NO_DATA__" >> $species_name"/Multimaping/Counts/"$species_name"_comparisons.tab"
  echo "No reads unicos todos: __REPLACE_NO_UNICO1__" >> $species_name"/Multimaping/Counts/"$species_name"_comparisons.tab"
  echo "No reads unicos, si con mapeo estandard: __REPLACE_NO_UNICO2__" >> $species_name"/Multimaping/Counts/"$species_name"_comparisons.tab"
  echo "No reads tradicional: __REPLACE_NO_Tradicional1__"  >> $species_name"/Multimaping/Counts/"$species_name"_comparisons.tab"
  echo "No reads tradicional pero si unicos: __REPLACE_NO_Tradicional2__"  >> $species_name"/Multimaping/Counts/"$species_name"_comparisons.tab"

  while [ $count -le $recorrer ]
  do
    imp_ID=$(sed -n $count"p" $imparalog_file | awk -F "\t" '{print $1}')
    genes=$(sed -n $count"p" $imparalog_file | awk -F "\t" '{print $2}' | tr -d "," | tr " " "\n")

    echo $imp_ID": "$genes
    for gen in $genes
    do
      original=$(grep $gen $species_name"/Work_Tables/"$species_name"_all_counts.tab" | tr "\t" "\n" | sed "1d" | awk '{sum+=$1;} END {print sum;}' | awk '{if ($1 >= 1) print "Data"; else print "Mising"}')
      multiple=$(grep $gen $species_name"/Multimaping/Counts/"$species_name"_all_counts.tab" | tr "\t" "\n" | sed "1d" | awk '{sum+=$1;} END {print sum;}' | awk '{if ($1 >= 1) print "Data"; else print "Mising"}')

      if [ "$original" == "Mising" ]
      then
        one_mapping_zero_count=$(($one_mapping_zero_count + 1))
      fi

      if [ "$multiple" == "Mising" ]
      then
        multimapping_mapping_zero_count=$(($multimapping_mapping_zero_count + 1))
      fi

      if [ "$original" == "Data" ] && [ "$multiple" == "Data" ]
      then
        type=$(echo "Consistent data")
      elif [ "$original" == "Mising" ] && [ "$multiple" == "Mising" ]
      then
        type=$(echo "No Mappings")
        all_missing=$(($all_missing + 1))
      elif [ "$original" == "Mising" ]
      then
        type=$(echo "No unique mappings")
        no_multi_mappings=$(($no_multi_mappings + 1))
      elif [ "$multiple" == "Mising" ]
      then
        type=$(echo "What?!!")
        weirdness=$(($weirdness + 1))
      fi
    done
    count=$(($count + 1))

    echo $imp_ID";"$gen";"$type | tr ";" "\t" >> $species_name"/Multimaping/Counts/"$species_name"_comparisons.tab"
  done

  sed -i "s/__REPLACE_NO_DATA__/$all_missing/" $species_name"/Multimaping/Counts/"$species_name"_comparisons.tab"
  sed -i "s/__REPLACE_NO_UNICO1__/$one_mapping_zero_count/" $species_name"/Multimaping/Counts/"$species_name"_comparisons.tab"
  sed -i "s/__REPLACE_NO_UNICO2__/$no_multi_mappings/" $species_name"/Multimaping/Counts/"$species_name"_comparisons.tab"
  sed -i "s/__REPLACE_NO_Tradicional1__/$multimapping_mapping_zero_count/" $species_name"/Multimaping/Counts/"$species_name"_comparisons.tab"
  sed -i "s/__REPLACE_NO_Tradicional2__/$weirdness/" $species_name"/Multimaping/Counts/"$species_name"_comparisons.tab"
}

get_TMM_data () {
  imparalog=$(sed -n $count"p"  "New_Sup_Tables/Temp/imparalog_genes_for_tmm.tmp" | awk -F ";" '{print $1}')
  gen=$(sed -n $count"p"  "New_Sup_Tables/Temp/imparalog_genes_for_tmm.tmp" | awk -F ";" '{print $2}')
  paper=$(sed -n $count"p"  "New_Sup_Tables/Temp/imparalog_genes_for_tmm.tmp" | awk -F ";" '{print $3}')
  study=$(sed -n $count"p"  "New_Sup_Tables/Temp/imparalog_genes_for_tmm.tmp" | awk -F ";" '{print $4}')
  stage=$(sed -n $count"p"  "New_Sup_Tables/Temp/imparalog_genes_for_tmm.tmp" | awk -F ";" '{print $5}')
  run=$(sed -n $count"p"  "New_Sup_Tables/Temp/imparalog_genes_for_tmm.tmp" | awk -F ";" '{print $6}')
  data=$(sed -n $count"p"  "New_Sup_Tables/Temp/imparalog_genes_for_tmm.tmp" | awk -F ";" '{print $7}')

  run_column=$(head -n 1 "Analisis_R/"$data"/"$data"_norm_counts.tab" | tr "\t" "\n" | grep -n $run | awk -F ":" '{print $1 + 1}' )
  tmm=$(grep $gen "Analisis_R/"$data"/"$data"_norm_counts.tab" | awk -v N=$run_column -F "\t" '{print $N}' | awk '{printf "%.2f\n", $1}')

  echo "Debug imparalog: "$imparalog
  echo "Debug gen: "$gen
  echo "Debug paper: "$paper
  echo "Debug study: "$study
  echo "Debug stage: "$stage
  echo "Debug run: "$run
  echo "Debug run_column: "$run_column
  echo "Debug tmm: "$tmm
  echo ""
}

################################################################################
################################################################################
################################################################################

if [ "$Schisto_TRIM" == "TRUE" ]
then
  species_name=Schistosoma
  rm -r $species_name"/raw_reads"

  mkdir $species_name"/raw_reads"
  sra_files=$(ls $species_name"/"*".sra" | sed 's/.sra//' | sed "s/.*\///")
  for sra in $sra_files
  do
    echo $sra
    fastq-dump --defline-seq '@$sn[_$rn]/$ri' --split-files $species_name"/"$sra".sra" -O $species_name"/raw_reads"
    echo "Comprimiendo"
    gzip $species_name"/raw_reads/"$sra"_1.fastq"
    gzip $species_name"/raw_reads/"$sra"_2.fastq"
  done
  run_trim
fi

################################################################################
################################################################################
################################################################################

if [ "$Schisto_STAR" == "TRUE" ]
then
  species_name=Schistosoma
  genome_fasta=$schisto_genome_fasta
  gff_file=$schisto_gff_file
  star_aln
fi

################################################################################
################################################################################
################################################################################

if [ "$Schisto_COUNT" == "TRUE" ]
then
  species_name=Schistosoma
  mrna=$schisto_trans
  gff_shenanigas
  counting
fi

################################################################################
################################################################################
################################################################################

if [ "$Schisto_TABLE_MAKER" == "TRUE" ]
then
  species_name=Schistosoma
  table_make
fi

################################################################################
################################################################################
################################################################################

if [ "$Schisto_Protasio_TRIM" == "TRUE" ]
then
  species_name=Schistosoma_Protasio2012
  rm -r $species_name"/raw_reads"

  mkdir $species_name"/raw_reads"
  sra_files=$(ls $species_name"/"*".sra" | sed 's/.sra//' | sed "s/.*\///")
  for sra in $sra_files
  do
    echo $sra
    fastq-dump --defline-seq '@$sn[_$rn]/$ri' --split-files $species_name"/"$sra".sra" -O $species_name"/raw_reads"
    echo "Comprimiendo"
    gzip $species_name"/raw_reads/"$sra"_1.fastq"
    gzip $species_name"/raw_reads/"$sra"_2.fastq"
  done
  run_trim
fi

################################################################################
################################################################################
################################################################################

if [ "$Schisto_Protasio_STAR" == "TRUE" ]
then
  species_name=Schistosoma_Protasio2012
  genome_fasta=$schisto_genome_fasta
  gff_file=/home/amanda/PROYECTS/Mauri_Imparalogos/Schistosoma/work_Schistosoma.gtf

  star_aln
fi

################################################################################
################################################################################
################################################################################

if [ "$Schisto_Protasio_COUNT" == "TRUE" ]
then
  species_name=Schistosoma_Protasio2012
  mrna=$schisto_trans
  gff_shenanigas
  counting
fi

################################################################################
################################################################################
################################################################################

if [ "$Schisto_Protasio_TABLE_MAKER" == "TRUE" ]
then
  species_name=Schistosoma_Protasio2012
  table_make
fi

################################################################################
################################################################################
################################################################################

if [ "$Schisto_Sanger_TRIM" == "TRUE" ]
then
  species_name=Smansoni_SangerInstitute
  rm -r $species_name"/raw_reads"

  mkdir $species_name"/raw_reads"
  sra_files=$(ls $species_name"/"*".sra" | sed 's/.sra//' | sed "s/.*\///")
  for sra in $sra_files
  do
    echo $sra
    fastq-dump --defline-seq '@$sn[_$rn]/$ri' --split-files $species_name"/"$sra".sra" -O $species_name"/raw_reads"
    echo "Comprimiendo"
    gzip $species_name"/raw_reads/"$sra"_1.fastq"
    gzip $species_name"/raw_reads/"$sra"_2.fastq"
  done
  run_trim
fi

################################################################################
################################################################################
################################################################################

if [ "$Schisto_Sanger_STAR" == "TRUE" ]
then
  species_name=Smansoni_SangerInstitute
  genome_fasta=$schisto_genome_fasta
  gff_file=/home/amanda/PROYECTS/Mauri_Imparalogos/Schistosoma/work_Schistosoma.gtf
  star_aln
fi

################################################################################
################################################################################
################################################################################

if [ "$Schisto_Sanger_COUNT" == "TRUE" ]
then
  species_name=Smansoni_SangerInstitute
  mrna=$schisto_trans
  gff_shenanigas
  counting
fi

################################################################################
################################################################################
################################################################################

if [ "$Schisto_Sanger_TABLE_MAKER" == "TRUE" ]
then
  species_name=Smansoni_SangerInstitute
  table_make
fi

################################################################################
################################################################################
################################################################################

if [ "$Smansoni_Genome_Location" == "TRUE" ]
then
  species_name=Schistosoma
  imparalog_file=Smansoni_imparalogs.ids
  gff_file=$schisto_gff_file
  new_annotation=$schisto_trans
  genome_context
fi

if [ "$Smansoni_Destance_groups" == "TRUE" ]
then
  species_name=Schistosoma
  imparalog_file=Smansoni_imparalogs.ids
  gff_file=$schisto_gff_file
  correlations_by_proximity
fi

################################################################################
################################################################################
################################################################################

if [ "$Dif_Exp_Table" == "TRUE" ]
then
  rm $expresion_file"/"*"_expresion_heatmap.txt"
  expresion_file=/home/amanda/PROYECTS/Mauri_Imparalogos/Analisis_R
  datasets=$(ls -d $expresion_file"/"* | sed "s/.*\///g")
  imparalog_file=Smansoni_imparalogs.ids

  for data in $datasets
  do
    echo $data
    comparisons=$(ls $expresion_file"/"$data"/Results_"*"_2023-02-23.tab"  | sed "s/.*\///g" | sed "s/Results_//" | sed "s/_2023-02-23.tab//")
    encabezado=$(echo "ImparalogID GenID")
    echo $encabezado
    for comp in $comparisons
    do
      first_group=$(echo $comp | tr "_" "\n" | head -n 1)
      second_group=$(echo $comp | tr "_" "\n" | tail -n 1)

      print_first=$(grep -w $first_group $expresion_file"/"$data"/Dif_exp_Keystone.txt" | awk -F ";" '{print $2}')
      print_second=$(grep -w $second_group $expresion_file"/"$data"/Dif_exp_Keystone.txt" | awk -F ";" '{print $2}')

      encabezado=$(echo $encabezado" "$print_first"_vs_"$print_second)
    done

    count=1
    recorrer=$(grep -c . $imparalog_file)
    echo $encabezado | tr " " "\t" >> $expresion_file"/"$data"_expresion_heatmap.txt"

    while [ $count -le $recorrer ]
    do
      echo $count"/"$recorrer

      imp_ID=$(sed -n $count"p" $imparalog_file | awk -F "\t" '{print $1}')
      genes=$(sed -n $count"p" $imparalog_file | awk -F "\t" '{print $2}' | tr -d ",")

      for gen in $genes
      do
        reportar=$(echo $imp_ID" "$gen)

        for comp in $comparisons
        do
          check_data=$(grep -c -w $gen $expresion_file"/"$data"/Results_"$comp"_2023-02-23.tab")

          if [ $check_data -eq 0 ]
          then
            expresion_change=NA
          else
            FDR=$(grep -w $gen $expresion_file"/"$data"/Results_"$comp"_2023-02-23.tab" | awk -F "\t" '{if ($6<=0.05) print "Sig"; else print "NoSig"}')
            Neg_Log_Fold_Change=$(grep -w $gen $expresion_file"/"$data"/Results_"$comp"_2023-02-23.tab" | awk -F "\t" '{if ($2>=0) print "Pos"; else print "Neg"}')

            if [ "$FDR" == "NoSig" ]
            then
              expresion_change=NoSig
            elif [ "$Neg_Log_Fold_Change" == "Pos" ]
            then
              expresion_change=$(grep -w $gen $expresion_file"/"$data"/Results_"$comp"_2023-02-23.tab" | awk -F "\t" '{if ($2>=2) print "Up"; else print "Stable"}')
            else
              expresion_change=$(grep -w $gen $expresion_file"/"$data"/Results_"$comp"_2023-02-23.tab" | awk -F "\t" '{if ($2<=-2) print "Down"; else print "Stable"}')
            fi
          fi
          reportar=$(echo $reportar" "$expresion_change)
        done
        echo $reportar | tr " " "\t" >> $expresion_file"/"$data"_expresion_heatmap.txt"
      done
      count=$(($count + 1))
    done
  done
fi

################################################################################
################################################################################
################################################################################

if [ "$Final_Summary" == "TRUE" ]
then
  wanted_data=$(echo "Schistosoma Schistosoma_Protasio2012 Smansoni_SangerInstitute")

  rm -r Paper_Sup_Material
  mkdir Paper_Sup_Material

  for data in $wanted_data
  do
    echo $data
    # seqkit stats -a -T $data"/raw_reads/"*".fastq.gz" >> "Paper_Sup_Material/"$data"_raw_reads_stats.txt"
    cp $data"/Work_Tables/"$data"_all_counts.tab" "Paper_Sup_Material/"$data"_all_counts.tab"

    read_files=$(ls $data"/STAR_aln/"*"_alnLog.final.out" | sed "s/.*\///" | sed "s/_alnLog.final.out//")
    count=6

    echo "Item "$read_files | tr " " "\t" >> "Paper_Sup_Material/"$data"_Alignment_Statistics.txt"
    while [ $count -le 36 ]
    do
      echo "inner while: "$count
      registro=$(sed -n $count"p" $data"/STAR_aln/"*"_alnLog.final.out" | awk -F "|" '{print $1}' | tr -s " " | sed "s/^ //" | sort -u)
      for read in $read_files
      do
        current_value=$(sed -n $count"p" $data"/STAR_aln/"$read"_alnLog.final.out" | awk -F "|" '{print $2}' | tr -d "\t"  )
        registro=$(echo $registro";"$current_value)
      done
      echo $registro | tr ";" "\t" >> "Paper_Sup_Material/"$data"_Alignment_Statistics.txt"

      count=$(($count + 1))
    done
  done
fi

################################################################################
################################################################################
################################################################################

if [ "$confirmacion_recorrida" == "TRUE" ]
then
  N_nuevos_genes=$(grep -c . verificacion_viejo_nuevo/nuevos_genes.work)
  N_viejos_genes=$(grep -c . verificacion_viejo_nuevo/viejos_genes.work)

  compartidos_viejo=$(grep -c -f verificacion_viejo_nuevo/nuevos_genes.work verificacion_viejo_nuevo/viejos_genes.work)
  compartidos_nuevo=$(grep -c -f verificacion_viejo_nuevo/viejos_genes.work verificacion_viejo_nuevo/nuevos_genes.work)

  grep -v -f verificacion_viejo_nuevo/viejos_genes.work verificacion_viejo_nuevo/nuevos_genes.work > verificacion_viejo_nuevo/Nuevas_adquisiciones.ids
  grep -v -f verificacion_viejo_nuevo/nuevos_genes.work verificacion_viejo_nuevo/viejos_genes.work > verificacion_viejo_nuevo/Perdidos.ids

  echo Perdidos_anotacion
  Perdidos_anotacion=$(awk -F "\t" '{if ($3=="gene") print}' $schisto_gff_file | grep -c -f verificacion_viejo_nuevo/Perdidos.ids)
  echo nuevos_anotacion
  nuevos_anotacion=$(awk -F "\t" '{if ($3=="gene") print}' $schisto_gff_file | grep -c -f verificacion_viejo_nuevo/nuevos_genes.work)
  echo viejos_anotacion
  viejos_anotacion=$(awk -F "\t" '{if ($3=="gene") print}' $schisto_gff_file | grep -c -f verificacion_viejo_nuevo/viejos_genes.work)

  echo "Total Nuevos: "$N_nuevos_genes > verificacion_viejo_nuevo/Resumen.txt
  echo "Total Viejos: "$N_viejos_genes >> verificacion_viejo_nuevo/Resumen.txt
  echo "compartidos N/V: "$compartidos_nuevo"/"$compartidos_nuevo >> verificacion_viejo_nuevo/Resumen.txt
  echo ""  >> verificacion_viejo_nuevo/Resumen.txt
  echo "Nuevos anotados: "$nuevos_anotacion >> verificacion_viejo_nuevo/Resumen.txt
  echo "Viejos anotados: "$viejos_anotacion >> verificacion_viejo_nuevo/Resumen.txt
  echo "Perdidos anotados: "$Perdidos_anotacion >> verificacion_viejo_nuevo/Resumen.txt
fi

################################################################################
################################################################################
################################################################################

if [ "$INTERPRO_PREPARE" == "TRUE" ]
then
  rm -r Interpro_Analisys/Input_files
  mkdir Interpro_Analisys/Input_files

  count=2
  recorrer=$(grep -c . $run_interpro_imparalogs_main_result"/inparalogs_group_composition.tsv")

  while [ $count -le $recorrer ]
  do
    echo $count"/"$recorrer

    genes=$(sed -n $count"p" $run_interpro_imparalogs_main_result"/inparalogs_group_composition_final.tsv" | awk -F "\t" '{print $1}' | tr -d " " | tr "," "\n")
    species_raw=$(sed -n $count"p" $run_interpro_imparalogs_main_result"/inparalogs_group_composition_final.tsv" | awk -F "\t" '{print $3}' | tr -d " ")
    species_work=$(grep -F -w $species_raw Species_keystone.ids | awk -F "\t" '{print $1}' )
    protein_file=$(ls $run_interpro_sequence_dir"/"$species_work"."*)

    for gen in $genes
    do
      internal_id=$(grep -w -F $gen $run_interpro_imparalogs_main_result"/gene_code_correspondance.tsv" | awk -F "\t" '{print $2}')
      control_extracton=$(seqkit grep -p $internal_id $protein_file | grep -c .)

      seqkit grep -p $internal_id $protein_file | sed "s/$internal_id/$gen/" >> "Interpro_Analisys/Input_files/"$species_work"_work.fasta"
    done
    count=$(($count + 1))
  done
fi

################################################################################
################################################################################
################################################################################

if [ "$INTERPRO_RUN" == "TRUE" ]
then
  rm -r "Interpro_Analisys/Interpro_Results"
  mkdir "Interpro_Analisys/Interpro_Results"

  especies=$(ls "Interpro_Analisys/Input_files/"*"_work.fasta" | sed "s/.*\///" | sed "s/_work.fasta//")
  for spe in $especies
  do
    echo "Running: "$spe
    $interprot_path"/interproscan.sh" -f TSV  -b "Interpro_Analisys/Interpro_Results/"$spe"_interpro" -dra -i "Interpro_Analisys/Input_files/"$spe"_work.fasta"
  done
fi

################################################################################
################################################################################
################################################################################

if [ "$INTERPRO_DOMAINS_IMP" == "TRUE" ]
then
  echo "Species Imp_ID Interpro_Domains" | tr " " "\t" > "Interpro_Analisys/Imparalog_domain_data.txt"

  count=2
  recorrer=$(grep -c . $run_interpro_imparalogs_main_result"/inparalogs_group_composition_final.tsv")
  while [ $count -le $recorrer ]
  do
    echo "Monitoreo: "$count"/"$recorrer

    genes=$(sed -n $count"p" $run_interpro_imparalogs_main_result"/inparalogs_group_composition_final.tsv" | awk -F "\t" '{print $1}' | tr -d " " | tr "," "\n")
    ImpID=$(sed -n $count"p" $run_interpro_imparalogs_main_result"/inparalogs_group_composition_final.tsv" | awk -F "\t" '{print $2}')
    species_raw=$(sed -n $count"p" $run_interpro_imparalogs_main_result"/inparalogs_group_composition_final.tsv" | awk -F "\t" '{print $3}' | tr -d " ")
    species_work=$(grep -F -w $species_raw "Species_keystone.ids" | awk -F "\t" '{print $1}' )

    registro_dominios=NA
    abemus_data=FALSE
    for gen in $genes
    do
      domains=$(awk -F "\t" -v GenID=$gen '{if ($1==GenID) print $12}' "Interpro_Analisys/Interpro_Results/"$species_work"_interpro.tsv" | grep IPR | sort -u)
      check_results=$(echo $domains | tr " " "\n" | grep -c .)

      if [ $check_results -gt 0 ]
      then
        abemus_data=TRUE
        registro_dominios=$(echo $registro_dominios" "$domains)
      fi
    done

    if [ "$abemus_data" == "TRUE" ]
    then
      registro_dominios=$(echo $registro_dominios | tr " " "\n" | grep -v -w NA | sort -u | tr "\n" ";" | sed "s/;$//")
    fi
    echo $species_work" "$ImpID" "$registro_dominios | tr " " "\t" >> "Interpro_Analisys/Imparalog_domain_data.txt"

    count=$(($count +1))
  done

  awk -F "\t"  '{print $12": "$13}' "Interpro_Analisys/Interpro_Results/"*"_interpro.tsv" | grep IPR | sort -u > "Interpro_Analisys/Domain_details.txt"
fi

################################################################################
################################################################################
################################################################################

if [ "$INTERPRO_GO_TERMS" == "TRUE" ]
then
  rm -r "Interpro_Analisys/GOTERM_Reference"
  mkdir "Interpro_Analisys/GOTERM_Reference"
  mkdir "Interpro_Analisys/GOTERM_Reference/Temp"

  echo "Especie ImparalogID GO_Terms" | tr " " "\t" >> "Interpro_Analisys/GOTERM_Reference/Imparalog_goterms.txt"

  especies=$(awk -F "\t" '{print $1}' "Species_keystone.ids")
  for spe in $especies
  do
    search_species=$(grep -w -F $spe "Species_keystone.ids" | awk -F "\t" '{print $2}' | awk -F "." '{print $2}')
    imparalogIDs=$(grep -w -F $search_species $run_interpro_imparalogs_main_result"/inparalogs_group_composition_final.tsv" | awk -F "\t" '{print $2}')

    grep -w -F $search_species $run_interpro_imparalogs_main_result"/inparalogs_group_composition_final.tsv" >> "Interpro_Analisys/GOTERM_Reference/Temp/"$spe"_imparalog_data.tmp"
    grep -w -F $search_species "Interpro_Analisys/GOTERM_Reference/Temp/"$spe"_imparalog_data.tmp" | awk -F "\t" '{print $1}' | tr -d " " | tr "," "\n" >> "Interpro_Analisys/GOTERM_Reference/Temp/"$spe"_all_original_ids.tmp"
    grep -w -F -f "Interpro_Analisys/GOTERM_Reference/Temp/"$spe"_all_original_ids.tmp" $run_interpro_imparalogs_main_result"/gene_code_correspondance.tsv" | awk -F "\t" '{print $2}' >> "Interpro_Analisys/GOTERM_Reference/Temp/"$spe"_all_internal_ids.tmp"
    grep -w -F -f "Interpro_Analisys/GOTERM_Reference/Temp/"$spe"_all_internal_ids.tmp" $eggnog_annotation | awk -F "\t" '{print $1"\t"$13}' >> "Interpro_Analisys/GOTERM_Reference/Temp/"$spe"_eggnog.tmp"

    for imp in $imparalogIDs
    do
      echo $spe" -- "$imp
      genes=$(grep -w -F $imp "Interpro_Analisys/GOTERM_Reference/Temp/"$spe"_imparalog_data.tmp" | awk -F "\t" '{print $1}' | tr -d " " | tr "," "\n")
      go_data=FALSE
      registro_go=NA
      for gen in $genes
      do
        internal_gen=$(grep -w -F $gen  $run_interpro_imparalogs_main_result"/gene_code_correspondance.tsv" | awk -F "\t" '{print $2}')
        check_go=$(awk -F "\t" -v gen=$internal_gen '{if ($1==gen) print $2}' "Interpro_Analisys/GOTERM_Reference/Temp/"$spe"_eggnog.tmp" | grep -c GO)

        if [ $check_go -gt 0 ]
        then
          go_data=TRUE
          gene_gos=$(awk -F "\t" -v gen=$internal_gen '{if ($1==gen) print $2}' "Interpro_Analisys/GOTERM_Reference/Temp/"$spe"_eggnog.tmp")
          registro_go=$(echo $registro_go","$gene_gos)
        fi
      done

      if [ "$go_data" == "TRUE" ]
      then
        registro_go=$(echo $registro_go | tr "," "\n" | grep -v NA | sort -u | tr "\n" "," | sed "s/,$//")
      fi

      echo $spe" "$imp" "$registro_go | tr " " "\t" >> "Interpro_Analisys/GOTERM_Reference/Imparalog_goterms.txt"
    done
  done
fi


################################################################################
################################################################################
################################################################################

if [ "$INTERPRO_DOMAINS_PRE_HEATMAP" == "TRUE" ]
then
  rm -r "Interpro_Analisys/Heatmap_Files"
  mkdir "Interpro_Analisys/Heatmap_Files"
  mkdir "Interpro_Analisys/Heatmap_Files/Temp"

  # Go_term_list=$(cat Interesting_GO_Terms.in)
  Go_term_list=$(echo "GO:0016787 GO:0008233 GO:0008289 GO:0016491 GO:0006952 GO:0043167 GO:0022857 GO:0016740")
  species_order=$(echo "hymenolepis_diminuta hymenolepis_microstoma mesocestoides_corti echinococcus_granulosus echinococcus_multilocularis echinococcus_canadensis taenia_solium taenia_asiatica taenia_saginata clonorchis_sinensis opistorchis_viverrini fasciola_hepatica trichobilharzia_regenti schistosoma_curassoni schistosoma_haematobium schistosoma_japonicum schistosoma_mansoni schistosoma_margrebowiei schistosoma_mattheei macrostomum_lignano schmidtea_mediterranea")
  domain_types=$(echo "Active_site Binding_site Conserved_site Domain Family Homologous_superfamily PTM Repeat")
  try_cut_offs=$(echo "0 1")

  for cut_off in $try_cut_offs
  do
    for go in $Go_term_list
    do
      print_go=$(echo $go | tr ":" "_")
      especies=$(grep -w -F $go "Interpro_Analisys/GOTERM_Reference/Imparalog_goterms.txt" | awk -F "\t" '{print $1}' | sort -u)
      total_especies=$(echo $especies | tr " " "\n" | grep -c .)
      grep -w -F $go "Interpro_Analisys/GOTERM_Reference/Imparalog_goterms.txt" | awk -F "\t" '{print $2}' > "Interpro_Analisys/Heatmap_Files/Temp/imparalog.tmp"
      all_domains=$(grep -w -F -f "Interpro_Analisys/Heatmap_Files/Temp/imparalog.tmp" "Interpro_Analisys/Imparalog_domain_data.txt" | awk -F "\t" '{print $3}' | tr ";" "\n" | sort -u | grep -v NA)
      for type in $domain_types
      do
        selected=NO_DATA
        for domain in $all_domains
        do
          check_good_domain=$(grep -F 'interpro id="'$domain'"' $interpro_ID_info | grep -c 'type="'$type'"')
          if [ $check_good_domain -gt 0 ]
          then
            echo "Selecting "$go": "$domain
            n_species=$(grep -f "Interpro_Analisys/Heatmap_Files/Temp/imparalog.tmp" "Interpro_Analisys/Imparalog_domain_data.txt" | grep $domain | awk -F "\t" '{print $1}' | sort -u | grep -c .)
            if [ $n_species -gt $cut_off ]
            then
              selected=$(echo $selected";"$domain)
            fi
            echo $domain": "$n_species"/"$total_especies >> "Interpro_Analisys/Heatmap_Files/"$print_go"_total_"$type"_"$cut_off".txt"
          fi
        done

        selected=$(echo $selected | sed "s/NO_DATA//" | tr ";" "\n" | sort -u | grep .)
        fix_N_domains=$(echo $selected | tr " " "\n" | grep -c . | awk '{print $1+1}')

        echo "Species "$selected | tr " " "\t" >> "Interpro_Analisys/Heatmap_Files/"$print_go"_type_"$type"_heatmap_"$cut_off"_norm.in"
        echo "Species "$selected | tr " " "\t" >> "Interpro_Analisys/Heatmap_Files/"$print_go"_type_"$type"_heatmap_"$cut_off"_raw.in"

        for spe in $especies
        do
          spe_gene_identifier=$(grep -w $spe "Species_keystone.ids" | awk -F "\t" '{print $3}')
          go_term_gene_count=$(grep -w $go $eggnog_annotation | grep -c -P -w ^$spe_gene_identifier)
          go_term_imparalog_count=$(grep -w $go "Interpro_Analisys/GOTERM_Reference/Imparalog_goterms.txt" | grep -c -w $spe )

          print_result_norm=$spe
          print_result_raw=$spe
          for domain in $selected
          do
            echo "Collecting "$go": "$spe" --- "$domain
            check_presence_norm=$(grep -w -F $spe "Interpro_Analisys/Imparalog_domain_data.txt" | grep -w -f "Interpro_Analisys/Heatmap_Files/Temp/imparalog.tmp" | grep -w -c $domain | awk -v total=$go_term_imparalog_count '{printf "%.10f\n", $1/total}')
            check_presence_raw=$(grep -w -F $spe "Interpro_Analisys/Imparalog_domain_data.txt" | grep -w -f "Interpro_Analisys/Heatmap_Files/Temp/imparalog.tmp" | grep -w -c $domain)

            print_result_norm=$(echo $print_result_norm" "$check_presence_norm)
            print_result_raw=$(echo $print_result_raw" "$check_presence_raw)
          done

          echo $print_result_norm | tr " " "\t" >> "Interpro_Analisys/Heatmap_Files/"$print_go"_type_"$type"_heatmap_"$cut_off"_norm.in"
          echo $print_result_raw | tr " " "\t" >> "Interpro_Analisys/Heatmap_Files/"$print_go"_type_"$type"_heatmap_"$cut_off"_raw.in"
        done

        check_data=$(grep -c . "Interpro_Analisys/Heatmap_Files/"$print_go"_type_"$type"_heatmap_"$cut_off"_norm.in")
        if [ $check_data -gt 1 ]
        then
          for domain in $selected
          do
            replace=$(grep $domain "Interpro_Analisys/Domain_details.txt" | sed "s/: /.../" | tr "/" "_" | tr " " "_")
            echo "Renaming: "$domain" -- "$replace
            sed -i "s/$domain/$replace/" "Interpro_Analisys/Heatmap_Files/"$print_go"_type_"$type"_heatmap_"$cut_off"_norm.in"
            sed -i "s/$domain/$replace/" "Interpro_Analisys/Heatmap_Files/"$print_go"_type_"$type"_heatmap_"$cut_off"_raw.in"
          done

          check_sufficient_data=$(head -n 1 "Interpro_Analisys/Heatmap_Files/"$print_go"_type_"$type"_heatmap_"$cut_off"_norm.in" | tr "\t" "\n" | grep -c .)
          if [ $check_sufficient_data -gt 3 ]
          then
            head -n 1 "Interpro_Analisys/Heatmap_Files/"$print_go"_type_"$type"_heatmap_"$cut_off"_norm.in" >> "Interpro_Analisys/Heatmap_Files/"$print_go"_heatmap_"$type"_"$cut_off"_work_norm.in"
            head -n 1 "Interpro_Analisys/Heatmap_Files/"$print_go"_type_"$type"_heatmap_"$cut_off"_raw.in" >> "Interpro_Analisys/Heatmap_Files/"$print_go"_heatmap_"$type"_"$cut_off"_work_raw.in"
            for spe_ord in $species_order
            do
              echo "Sorting: "$print_go" -- "$spe_ord
              check_species=$(grep -w -c $spe_ord "Interpro_Analisys/Heatmap_Files/"$print_go"_type_"$type"_heatmap_"$cut_off"_norm.in" )
              if [ $check_species -gt 0 ]
              then
                grep -w $spe_ord "Interpro_Analisys/Heatmap_Files/"$print_go"_type_"$type"_heatmap_"$cut_off"_norm.in" >> "Interpro_Analisys/Heatmap_Files/"$print_go"_heatmap_"$type"_"$cut_off"_work_norm.in"
                grep -w $spe_ord "Interpro_Analisys/Heatmap_Files/"$print_go"_type_"$type"_heatmap_"$cut_off"_raw.in" >> "Interpro_Analisys/Heatmap_Files/"$print_go"_heatmap_"$type"_"$cut_off"_work_raw.in"
              else
                missing_domains=$(seq -s= $fix_N_domains | sed "s/=/\n=\n/g" | grep = | tr "=" "0" | tr "\n" " " | sed "s/$/\n/")
                echo $spe_ord" "$missing_domains | tr " " "\t"  >> "Interpro_Analisys/Heatmap_Files/"$print_go"_heatmap_"$type"_"$cut_off"_work_norm.in"
                echo $spe_ord" "$missing_domains | tr " " "\t"  >> "Interpro_Analisys/Heatmap_Files/"$print_go"_heatmap_"$type"_"$cut_off"_work_raw.in"
              fi
            done
          fi
        fi
      done
    done
  done
fi

################################################################################
################################################################################
################################################################################

if [ "$SUP_Tables" == "TRUE" ]
then
  rm -r "New_Sup_Tables"
  mkdir "New_Sup_Tables"
  mkdir "New_Sup_Tables/Temp"

  datasets=$(echo "Schistosoma Schistosoma_Protasio2012 Smansoni_SangerInstitute")

  # Table_Expresion_Correlation
  for data in $datasets
  do
    echo "Imparalog Gen1 Gen2 Dist_Classification Distance Used Pearson Speerman Interprot1 Interprot2" | tr " " "\t" >> "New_Sup_Tables/"$data"_Correlation_summary.txt"

    count=1
    recorrer_imparalogos=$(grep -c . "Analisis_R/"$data"/"$data"_Registro_uso.txt")

    while [ $count -le $recorrer_imparalogos ]
    do
      imparalog=$(sed -n $count"p" "Analisis_R/"$data"/"$data"_Registro_uso.txt" | awk -F "\t" '{print $2}' | tr -d '"')
      gen1=$(sed -n $count"p" "Analisis_R/"$data"/"$data"_Registro_uso.txt" | awk -F "\t" '{print $3}' | tr -d '"')
      gen2=$(sed -n $count"p" "Analisis_R/"$data"/"$data"_Registro_uso.txt" | awk -F "\t" '{print $4}' | tr -d '"')
      used=$(sed -n $count"p" "Analisis_R/"$data"/"$data"_Registro_uso.txt" | awk -F "\t" '{print $5}' | tr -d '"')

      distance=$(grep -w -F $imparalog "Analisis_R/"$data"/distance_corr.txt" | grep -w $gen1 | grep -w $gen2 | awk -F "\t" '{print $4"_oo_"$5}')
      correlation_info=$(grep -w -F $imparalog "Analisis_R/"$data"/"$data"_correlation_final.tab" | grep -w $gen1 | grep -w $gen2 | awk -F "\t" '{print $5"_oo_"$6}')

      interpro_info1=$(grep -w $gen1 "Interpro_Analisys/Interpro_Results/schistosoma_mansoni_interpro.tsv" | awk -F "\t" '{print $12" ("$13")"}' | sed "s/- (-)/NA/" | grep IPR | tr "\n" ";" | sed "s/;$/\n/")
      check_inter1=$(echo $interpro_info1 | grep -c IPR)
      if [ $check_inter1 -eq 0 ]
      then
        interpro_info1=NA
      fi

      interpro_info2=$(grep -w $gen2 "Interpro_Analisys/Interpro_Results/schistosoma_mansoni_interpro.tsv" | awk -F "\t" '{print $12" ("$13")"}' | sed "s/- (-)/NA/" | grep IPR | tr "\n" ";" | sed "s/;$/\n/")
      check_inter2=$(echo $interpro_info2 | grep -c IPR)
      if [ $check_inter2 -eq 0 ]
      then
        interpro_info2=NA
      fi
      echo "Imparalog: "$imparalog
      echo "Gen1: "$gen1
      echo "Gen2: "$gen2
      echo "Distance: "$distance
      echo "used: "$used
      echo "Correlation_info: "$correlation_info
      echo "interpro_info1: "$interpro_info1
      echo "interpro_info2: "$interpro_info2
      echo "########################################"
      echo ""

      echo $imparalog"_oo_"$gen1"_oo_"$gen2"_oo_"$distance"_oo_"$used"_oo_"$correlation_info"_oo_"$interpro_info1"_oo_"$interpro_info2 | sed "s/_oo_/\t/g" >> "New_Sup_Tables/"$data"_Correlation_summary.txt"

      count=$(($count + 1))
    done
  done

  # Preparar Tabla 3
  for data in $datasets
  do
    echo "Preparando: "$data
    if [ $data == "Schistosoma" ]
    then
      paper=$(echo "Wangwiwatsin et al. 2020")
      study=ERP113121
    elif [ "$data" == "Schistosoma_Protasio2012" ]
    then
      paper=$(echo "Protasio et al. 2012")
      study=ERP000427
    elif [ "$data" == "Smansoni_SangerInstitute" ]
    then
      paper=$(echo "Wellcome Sanger Institute" )
      study=ERP017466
    fi

    grupos=$(sed 1d "Analisis_R/"$data"/"$data"_grupos.tab" | awk -F "\t" '{print $2}' | sort -u)

    count=1
    recorrer_raw_imparalogos=$(grep -c . Smansoni_imparalogs.ids)
    while [ $count -le $recorrer_raw_imparalogos ]
    do
      imparalog=$(sed -n $count"p" Smansoni_imparalogs.ids | awk -F "\t" '{print $1}')
      genes=$(sed -n $count"p" Smansoni_imparalogs.ids | awk -F "\t" '{print $2}' | tr -d " " | tr "," "\n" | awk -F "." '{print $1}')

      echo $data" -- "$imparalog": "$genes
      for gen in $genes
      do
        for gru in $grupos
        do
          corridas=$(grep -F -w $gru "Analisis_R/"$data"/"$data"_grupos.tab" | awk -F "\t" '{print $1}' | sort)
          print_gru=$(echo $gru | sed "s/^1_//" | sed "s/^2_//" | sed "s/^3_//" | sed "s/^4_//" | sed "s/^5_//" | sed "s/^6_//" | tr "_" " ")
          for run in $corridas
          do
            echo $imparalog";"$gen";"$paper";"$study";"$print_gru";"$run";"$data >> "New_Sup_Tables/Temp/imparalog_genes_for_tmm.tmp"
          done
        done
      done
      count=$(($count + 1))
    done
  done

  recorrer_nueva_tabla=$(grep -c . "New_Sup_Tables/Temp/imparalog_genes_for_tmm.tmp")
  echo "Code of Group of Inparalogs;Gene code;Source;Study;Developmental stage;Run;TMM per run;Median TMM among replicates" | tr ";" "\t" >> "New_Sup_Tables/New_SupTab_TMM_Data.txt"
  count=1
  get_TMM_data

  echo $imparalog";"$gen";"$paper";"$study";"$stage";"$run";"$tmm";REPLACE_ME!!" | tr ";" "\t" >> "New_Sup_Tables/New_SupTab_TMM_Data.txt"

  previous_stage=$stage
  registro_tmm=$tmm
  count=2

  while [ $count -le $recorrer_nueva_tabla ]
  do
    get_TMM_data

    if [ "$stage" == "$previous_stage" ]
    then
      registro_tmm=$(echo $registro_tmm" "$tmm)
      echo $imparalog";"$gen";"$paper";"$study";"$stage";"$run";"$tmm";REPLACE_ME!!" | tr ";" "\t" >> "New_Sup_Tables/New_SupTab_TMM_Data.txt"
    else
      previous_stage=$stage
      median=$(echo $registro_tmm | tr " " "\n" | grep . | datamash median 1)
      echo "Debug median: "$median
      echo ""
      sed -i "s/REPLACE_ME!!/$median/" "New_Sup_Tables/New_SupTab_TMM_Data.txt"
      registro_tmm=$tmm
      echo $imparalog";"$gen";"$paper";"$study";"$stage";"$run";"$tmm";REPLACE_ME!!" | tr ";" "\t" >> "New_Sup_Tables/New_SupTab_TMM_Data.txt"
    fi
    count=$(($count +1))
  done
  median=$(echo $registro_tmm | tr " " "\n" | grep . | datamash median 1)
  sed -i "s/REPLACE_ME!!/$median/" "New_Sup_Tables/New_SupTab_TMM_Data.txt"
fi


if [ $SUP_ANOTATION == "TRUE" ]
then
  count=2
  recorrer=$(grep -c . "Interpro_Analisys/Imparalog_domain_data.txt")

  echo "Especie Imparalogo N_Genes N_Domains N_Goterms GenesID Interpro_Signatures GOTerms" | tr " " "\t" > "New_Sup_Tables/Sup_Annotation_Results.txt"

  while [ $count -le $recorrer ]
  do
    echo $count" / "$recorrer
    species=$(sed -n $count"p" "Interpro_Analisys/Imparalog_domain_data.txt" | awk -F "\t" '{print $1}')
    imp=$(sed -n $count"p" "Interpro_Analisys/Imparalog_domain_data.txt" | awk -F "\t" '{print $2}')
    genes=$(grep -w -F $imp $run_interpro_imparalogs_main_result"/inparalogs_group_composition_final.tsv" | awk -F "\t" '{print $1}' | tr -d " " )
    N_genes=$(echo $genes | tr "," "\n" | grep -c .)

    pre_domains=$(sed -n $count"p" "Interpro_Analisys/Imparalog_domain_data.txt" | awk -F "\t" '{print $3}' | tr ";" " ")
    check_pre_domains=$(echo $pre_domains | tr " " "\n" | grep -w -c NA)

    N_domains=0
    if [ $check_pre_domains -gt 0 ]
    then
      registro_dominios=NA
    else
      registro_dominios=NA
      for dom in $pre_domains
      do
        N_domains=$(($N_domains + 1))
        print_domain=$(grep -w $dom "Interpro_Analisys/Domain_details.txt" | sed "s/$/)/" | sed "s/: / (/")
        registro_dominios=$(echo $registro_dominios";"$print_domain)
      done
      registro_dominios=$(echo $registro_dominios | sed "s/^NA;//" | sed "s/;/; /g")
    fi

    goterms=$(grep -w -F $imp "Interpro_Analisys/GOTERM_Reference/Imparalog_goterms.txt" | awk -F "\t" '{print $3}')
    if [ $goterms == "NA" ]
    then
      N_goterms=0
    fi

    N_goterms=$(echo $goterms | tr "," "\n" | grep -c -v NA)
    echo $species"_ooo_"$imp"_ooo_"$N_genes"_ooo_"$N_domains"_ooo_"$N_goterms"_ooo_"$genes"_ooo_"$registro_dominios"_ooo_"$goterms | sed "s/_ooo_/\t/g" >> "New_Sup_Tables/Sup_Annotation_Results.txt"

    count=$(($count + 1))
  done
fi


if [ "$SUP_ANOTATION_CONSULTA" == "TRUE" ]
then
  echo "Especie N_Imparalogos N_Genes N_Interpro N_GoTerms" | tr " " "\t" > "New_Sup_Tables/Sup_Annotation_consulta.txt"

  especies=$(sed 1d "New_Sup_Tables/Sup_Annotation_Results.txt" | awk -F "\t" '{print $1}' | sort -u)
  for spe in $especies
  do
    N_Imparalogos=$(grep -w $spe "New_Sup_Tables/Sup_Annotation_Results.txt" |  grep -c .)
    N_Genes=$(grep -w $spe "New_Sup_Tables/Sup_Annotation_Results.txt" | awk -F "\t" '{print $5}' | tr -d " " | sed "s/,/\n/" | grep -c .)
    N_Imparalogos_interpro=$(grep -w $spe "New_Sup_Tables/Sup_Annotation_Results.txt" | awk -F "\t" '{if ($7 != "NA") print $5}' | tr -d " " | sed "s/,/\n/" | grep -c .)
    N_Imparalogos_goterm=$(grep -w $spe "New_Sup_Tables/Sup_Annotation_Results.txt" | awk -F "\t" '{if ($8 != "NA") print $5}' | tr -d " " | sed "s/,/\n/" | grep -c .)

    echo $spe" "$N_Imparalogos" "$N_Genes" "$N_Imparalogos_interpro" "$N_Imparalogos_goterm | tr " " "\t" >> "New_Sup_Tables/Sup_Annotation_consulta.txt"
  done
fi
