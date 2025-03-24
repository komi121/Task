#!/bin/bash

# Задаем переменные для входных и выходных файлов.
REFERENCE_GENOME="/mnt/8Tb_new/arkom/BWA_reference/GRCh38.primary_assembly.genome.cleaned.fa"  # Путь к референсному геному
DATA_DIR="/mnt/8Tb_new/arkom/test_task-1"  # Папка с данными
OUTPUT_DIR="/mnt/8Tb_new/arkom/results_task-1"  # Папка для результатов
STAT_OUTPUT="${OUTPUT_DIR}/all_stat.txt"  # Файл для статистики
THREADS=4  # Количество потоков 

# Функция для обработки одного образца
process_sample() {
    local SAMPLE_DIR=$1
    local OUTPUT_DIR=$2
    local STAT_OUTPUT=$3
    local SAMPLE_NAME=$(basename "$SAMPLE_DIR")
    local R1=$(find "$SAMPLE_DIR" -name "*_read1.fastq.gz")  # Находим файл R1
    local R2=$(find "$SAMPLE_DIR" -name "*_read2.fastq.gz")  # Находим файл R2

    local OUTPUT_DIR_PER_SAMPLE="${OUTPUT_DIR}/${SAMPLE_NAME}" # получаем путь для сохранения результатов 

    echo "Анализ образца: $SAMPLE_NAME"
    mkdir -p "$OUTPUT_DIR_PER_SAMPLE" # Создаем папку

    # FastQC анализ
    echo "Анализ FastQC"
    fastqc -o "$OUTPUT_DIR_PER_SAMPLE" -t $THREADS "$R1" "$R2"


    # Выравнивание прочтений с использованием BWA
    echo "Шаг 1: Выравнивание прочтений BWA"
    bwa mem -t $THREADS "$REFERENCE_GENOME" "$R1" "$R2" > "${OUTPUT_DIR_PER_SAMPLE}/${SAMPLE_NAME}.sam"

    # Преобразование SAM в BAM
    echo "Шаг 2: Конвертируем SAM в BAM"
    samtools view -@ $THREADS -Sb "${OUTPUT_DIR_PER_SAMPLE}/${SAMPLE_NAME}.sam" > "${OUTPUT_DIR_PER_SAMPLE}/${SAMPLE_NAME}.bam"
    rm "${OUTPUT_DIR_PER_SAMPLE}/${SAMPLE_NAME}.sam" # Удаляем SAM файл

    # Сортировка BAM
    echo "Шаг 3: Сортировка BAM файла"
    samtools sort -@ $THREADS -o "${OUTPUT_DIR_PER_SAMPLE}/${SAMPLE_NAME}_sorted.bam" "${OUTPUT_DIR_PER_SAMPLE}/${SAMPLE_NAME}.bam"
    samtools index "${OUTPUT_DIR_PER_SAMPLE}/${SAMPLE_NAME}_sorted.bam" # Сразу добавляем индексы для bam файла

    samtools stats "${OUTPUT_DIR_PER_SAMPLE}/${SAMPLE_NAME}_sorted.bam" > "${OUTPUT_DIR_PER_SAMPLE}/${SAMPLE_NAME}_stats.txt" # Сразу генерируем статистику
    rm "${OUTPUT_DIR_PER_SAMPLE}/${SAMPLE_NAME}.bam" # Удаляем bam файл

    # Добавление групп чтения и удаление дубликатов
    echo "Шаг 4: Анализ дупликатов"

    picard AddOrReplaceReadGroups \
        --INPUT "${OUTPUT_DIR_PER_SAMPLE}/${SAMPLE_NAME}_sorted.bam" \
        --OUTPUT "${OUTPUT_DIR_PER_SAMPLE}/${SAMPLE_NAME}_rg.bam" \
        --RGID "$SAMPLE_NAME" \
        --RGLB lib1 \
        --RGPL ILLUMINA \
        --RGPU unit1 \
        --RGSM "$SAMPLE_NAME" # Добавляем нужную информацию в bam файл, иначе picard работать не будет

    picard MarkDuplicates \
        --INPUT "${OUTPUT_DIR_PER_SAMPLE}/${SAMPLE_NAME}_rg.bam" \
        --OUTPUT "${OUTPUT_DIR_PER_SAMPLE}/${SAMPLE_NAME}_dedup.bam" \
        --METRICS_FILE "${OUTPUT_DIR_PER_SAMPLE}/${SAMPLE_NAME}_dedup_metrics.txt" # Генерируем статистику по дупликатам

    picard CollectInsertSizeMetrics \
        I="${OUTPUT_DIR_PER_SAMPLE}/${SAMPLE_NAME}_rg.bam" \
        O="${OUTPUT_DIR_PER_SAMPLE}/${SAMPLE_NAME}_InsertSizeMetrics.txt" \
        H="${OUTPUT_DIR_PER_SAMPLE}/${SAMPLE_NAME}_InsertSizeHistogram.pdf" \
        M=0.5 # Генерируем статистику по вставкам 

    # Оценка основных метрик качества
    echo "Шаг 5: Cобираем статистику в один файл"


    # Общее количество ридов 
    local TOTAL_READS=$(grep 'SN	raw total sequences' "${OUTPUT_DIR_PER_SAMPLE}/${SAMPLE_NAME}_stats.txt" | awk '{print $5}')

    # Длина ридов 
    local READ_LENGTH=$(grep 'SN	average length' "${OUTPUT_DIR_PER_SAMPLE}/${SAMPLE_NAME}_stats.txt" | awk '{print $4}')

    # Процент перекрывающихся ридов
    local TOTAL_PAIRED=$(samtools view -c -f 2 "${OUTPUT_DIR_PER_SAMPLE}/${SAMPLE_NAME}_sorted.bam") 
    local OVERLAPPING=$(samtools view -f 2 "${OUTPUT_DIR_PER_SAMPLE}/${SAMPLE_NAME}_sorted.bam" | awk '$4 + length($10) > $8' | wc -l)
    local OVERLAP_PERCENT=$(echo "scale=2; $OVERLAPPING / $TOTAL_PAIRED * 100" | bc)

    # Процент адаптеров 
    unzip -q "${OUTPUT_DIR_PER_SAMPLE}/${SAMPLE_NAME}_read1_fastqc.zip" -d "${OUTPUT_DIR_PER_SAMPLE}"
    unzip -q "${OUTPUT_DIR_PER_SAMPLE}/${SAMPLE_NAME}_read2_fastqc.zip" -d "${OUTPUT_DIR_PER_SAMPLE}"

    local ADAPTER_PERCENT_R1=$(tail -n 2 "${OUTPUT_DIR_PER_SAMPLE}/${SAMPLE_NAME}_read1_fastqc/fastqc_data.txt" | head -n 1 | awk '{print $2}')
    local ADAPTER_PERCENT_R2=$(tail -n 2 "${OUTPUT_DIR_PER_SAMPLE}/${SAMPLE_NAME}_read2_fastqc/fastqc_data.txt" | head -n 1 | awk '{print $2}')

    rm -r "${OUTPUT_DIR_PER_SAMPLE}/${SAMPLE_NAME}_read1_fastqc"
    rm -r "${OUTPUT_DIR_PER_SAMPLE}/${SAMPLE_NAME}_read2_fastqc"

    # Процент дупликатов
    local DUPLICATE_PERCENT=$(grep '^lib1' "${OUTPUT_DIR_PER_SAMPLE}/${SAMPLE_NAME}_dedup_metrics.txt" | awk '{print $9}')

    # Процент неуникально сопоставленных ридов
    local TOTAL_MAPPED=$(grep 'SN	reads mapped:' "${OUTPUT_DIR_PER_SAMPLE}/${SAMPLE_NAME}_stats.txt" | awk '{print $4}')
    local NON_UNIQUE=$(grep 'SN	reads MQ0:' "${OUTPUT_DIR_PER_SAMPLE}/${SAMPLE_NAME}_stats.txt" | awk '{print $4}')
    local NON_UNIQUE_PERCENT=$(echo "scale=4; $NON_UNIQUE * 100 / $TOTAL_READS" | bc)

        # Средний, медианный и модальный размер вставки
    local INSERT_SIZE_MEAN=$(awk -F'\t' '/MEDIAN_INSERT_SIZE/ { getline; print $6 }' "${OUTPUT_DIR_PER_SAMPLE}/${SAMPLE_NAME}_InsertSizeMetrics.txt")
    local INSERT_SIZE_MEDIAN=$(awk -F'\t' '/MEDIAN_INSERT_SIZE/ { getline; print $1 }' "${OUTPUT_DIR_PER_SAMPLE}/${SAMPLE_NAME}_InsertSizeMetrics.txt")
    local INSERT_SIZE_MODE=$(awk -F'\t' '/MEDIAN_INSERT_SIZE/ { getline; print $2 }' "${OUTPUT_DIR_PER_SAMPLE}/${SAMPLE_NAME}_InsertSizeMetrics.txt")

    echo "$SAMPLE_NAME:$TOTAL_READS:$READ_LENGTH:$OVERLAP_PERCENT:$ADAPTER_PERCENT_R1:$ADAPTER_PERCENT_R2:$DUPLICATE_PERCENT:$NON_UNIQUE_PERCENT:$INSERT_SIZE_MEAN:$INSERT_SIZE_MEDIAN:$INSERT_SIZE_MODE" >> "$STAT_OUTPUT"
}

# Создание файла для статистики
mkdir -p "$OUTPUT_DIR"
touch "$STAT_OUTPUT"
echo "Sample:All_reads:Read_length:Overlapping_reads:Adapter_percent_R1:Adapter_percent_R2:Duplicate_percent:Non_unique:Insert_size_mean:Insert_size_median:Insert_size_mode" > "$STAT_OUTPUT"

# Основной цикл: проход по всем папкам с образцами
for SAMPLE_DIR in "${DATA_DIR}"/*; do
    if [[ -d "$SAMPLE_DIR" ]]; then
        process_sample "$SAMPLE_DIR" "$OUTPUT_DIR" "$STAT_OUTPUT"
    fi
done

# Генерация отчета MultiQC для всех образцов
echo "Step 7: Generating MultiQC report..."
multiqc "$OUTPUT_DIR" -o "${OUTPUT_DIR}/multiqc_report"


echo "All samples processed! MultiQC report saved in ${OUTPUT_DIR}/multiqc_report."