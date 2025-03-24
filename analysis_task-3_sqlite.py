import os
import glob
import sqlite3
import matplotlib.pyplot as plt
import seaborn as sns

def get_vcf_files(vcf_files_path):
    """
    Получает все VCF файлы в указанной директории.
    """
    print("Получаем список vcf файлов")
    vcf_files = glob.glob(os.path.join(vcf_files_path, "*.vcf"))
    for index in range(len(vcf_files)):
         vcf_files[index] = vcf_files[index].split("/")[-1]
    
    print(f"Обнаружено {len(vcf_files)} VCF файлов")
    return vcf_files


def check_folder(folder_name):
    """
    Создает папку, если она не существует.
    """
    if not os.path.exists(folder_name):
        os.makedirs(folder_name)  # Рекурсивно создает все папки в пути
        print(f"Папка '{folder_name}' успешно создана.")
    else:
        pass


def setup_database(path_output,data_base):
    """
    Создает базу данных SQLite и необходимые таблицы.
    """
    check_folder(path_output)
    conn = sqlite3.connect(f"{path_output}/{data_base}")
    cursor = conn.cursor()
    
    # Таблица для вариантов
    cursor.execute('''
    CREATE TABLE IF NOT EXISTS variants (
        variant_key TEXT PRIMARY KEY,
        chrom TEXT,
        pos TEXT,
        id TEXT,
        ref TEXT,
        alt TEXT,
        count INTEGER
    )
    ''')
    
    # Таблица для образцов
    cursor.execute('''
    CREATE TABLE IF NOT EXISTS samples (
        sample_name TEXT PRIMARY KEY,
        count_all INTEGER,
        count_unique INTEGER DEFAULT 0
    )
    ''')
    
    # Таблица для связи вариантов и образцов
    cursor.execute('''
    CREATE TABLE IF NOT EXISTS variant_samples (
        variant_key TEXT,
        sample_name TEXT,
        qual TEXT,
        filter_fields TEXT,
        info_fields TEXT,
        format_fields TEXT,
        sample_data TEXT,
        depth TEXT,
        ref_alt_depth TEXT,
        alt_depth TEXT,
        VAF TEXT,
        PRIMARY KEY (variant_key, sample_name),
        FOREIGN KEY (variant_key) REFERENCES variants (variant_key),
        FOREIGN KEY (sample_name) REFERENCES samples (sample_name)
    )
    ''')
    
    # Таблица для заголовков VCF
    cursor.execute('''
    CREATE TABLE IF NOT EXISTS vcf_headers (
        sample_name TEXT,
        header_line TEXT,
        FOREIGN KEY (sample_name) REFERENCES samples (sample_name)
    )
    ''')
    
    # Таблица для фильтров
    cursor.execute('''
    CREATE TABLE IF NOT EXISTS filters (
        sample_name TEXT,
        filter_name TEXT,
        count INTEGER,
        PRIMARY KEY (sample_name, filter_name),
        FOREIGN KEY (sample_name) REFERENCES samples (sample_name)
    )
    ''')
    
    conn.commit()
    conn.close()


def parse_vcf_files(vcf_files, vcf_files_path, path_output,data_base):
    """
    Парсит VCF файлы и сохраняет данные в базу SQLite.
    """
    print("Начинаем парсить VCF файлы")
    conn = sqlite3.connect(f"{path_output}/{data_base}")
    cursor = conn.cursor()
    count_sample = len(vcf_files)
    i = 1
    print(f"В обработку принято {len(vcf_files)} файлов")
    
    for vcf_file in vcf_files:
        print(f"Начинаем обработку {vcf_file}")
        sample_name = os.path.basename(vcf_file).split('.')[0]
        
        # Добавляем образец в базу
        cursor.execute('''
        INSERT OR REPLACE INTO samples (sample_name, count_all) VALUES (?, 0)
        ''', (sample_name,))
        
        with open(vcf_files_path + vcf_file, 'rt') as f:
            variant_count = 0
            for line in f:
                if line.startswith('#'):
                    cursor.execute('''
                    INSERT INTO vcf_headers (sample_name, header_line) VALUES (?, ?)
                    ''', (sample_name, line))
                    continue
                
                variant_count += 1
                fields = line.strip().split('\t')
                chrom = fields[0]
                pos = fields[1]
                id = fields[2]
                ref = fields[3]
                alt = fields[4]
                qual = fields[5]
                filter_fields = fields[6]
                info_fields = fields[7]
                format_fields = fields[8]
                sample_data = fields[9]
                
                # Получаем информацию о глубине покрытия
                depth = None
                ref_alt_depth = None
                alt_depth = None
                vaf = None
                
                format_items = format_fields.split(':')
                sample_items = sample_data.split(':')
                
                if 'DP' in format_items:
                    dp_idx = format_items.index('DP')
                    if dp_idx < len(sample_items):
                        depth = sample_items[dp_idx]
                
                if 'AD' in format_items:
                    ad_idx = format_items.index('AD')
                    if ad_idx < len(sample_items):
                        ref_alt_depth = sample_items[ad_idx].replace(",", ":")
                        alt_depth = ref_alt_depth.split(":")[1]
                        if int(alt_depth) == 0:
                            vaf = 0
                        else:
                            vaf = int(alt_depth) / int(depth)

                
                # Создаем ключ варианта
                variant_key = f"{chrom}_{pos}_{ref}_{alt}"
                
                # Проверяем, существует ли уже такой вариант
                cursor.execute('''
                SELECT count FROM variants WHERE variant_key = ?
                ''', (variant_key,))
                
                result = cursor.fetchone()
                
                if result:
                    # Обновляем счетчик варианта
                    count = result[0] + 1
                    cursor.execute('''
                    UPDATE variants SET count = ? WHERE variant_key = ?
                    ''', (count, variant_key))
                else:
                    # Добавляем новый вариант
                    cursor.execute('''
                    INSERT INTO variants (variant_key, chrom, pos, id, ref, alt, count)
                    VALUES (?, ?, ?, ?, ?, ?, 1)
                    ''', (variant_key, chrom, pos, id, ref, alt))

                # Добавляем связь варианта и образца
                cursor.execute('''
                INSERT INTO variant_samples 
                (variant_key, sample_name, qual, filter_fields, info_fields, format_fields, sample_data, depth, ref_alt_depth,alt_depth,VAF)
                VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?,?,?)
                ''', (variant_key, sample_name, qual, filter_fields, info_fields, format_fields, sample_data, depth, ref_alt_depth,alt_depth,vaf))
            
            # Обновляем количество вариантов для образца
            cursor.execute('''
            UPDATE samples SET count_all = ? WHERE sample_name = ?
            ''', (variant_count, sample_name))
        
        i += 1
        print(f"Обработку {vcf_file} закончил")
        print(f"Осталось {count_sample - i + 1} из {count_sample} образцов")
        
        # Сохраняем изменения после обработки каждого файла
        conn.commit()
    conn.close()

def generate_files(conn, folder_output, add_coment_in_output_file):
    """
    Генерирует файлы с уникальными вариантами для каждого образца.
    """
    check_folder(folder_output)
    cursor = conn.cursor()
    
    # Получаем все образцы
    cursor.execute('''SELECT sample_name FROM samples''')
    samples = cursor.fetchall()
    
    for sample in samples:
        sample_name = sample[0]
        
        with open(f"{folder_output}/{sample_name}_unique.vcf", 'w') as f:
            if add_coment_in_output_file:
                # Получаем все заголовки для образца
                cursor.execute('''
                SELECT header_line FROM vcf_headers WHERE sample_name = ?
                ''', (sample_name,))
                
                headers = cursor.fetchall()
                for header in headers:
                    f.write(header[0])
            else:
                # Получаем только последний заголовок
                cursor.execute('''
                SELECT header_line FROM vcf_headers WHERE sample_name = ?
                ORDER BY rowid DESC LIMIT 1
                ''', (sample_name,))
                
                header = cursor.fetchone()
                if header:
                    f.write(header[0])


def write_unique_variants( folder_output, add_coment_in_output_file,data_base ):
    """
    Записывает уникальные варианты в файлы.
    """
    print("Записываем уникальные snp для каждого варианта")
    conn = sqlite3.connect(f"{folder_output}/{data_base}")
    cursor = conn.cursor()
    
    # Генерация файлов с уникальными вариантами
    path_output = f"{folder_output}/unique_variants"
    generate_files(conn, path_output, add_coment_in_output_file)
    
    # Получаем уникальные варианты (count = 1)
    cursor.execute('''
    SELECT v.variant_key, v.chrom, v.pos, v.id, v.ref, v.alt, 
           vs.sample_name, vs.qual, vs.filter_fields, vs.info_fields, 
           vs.format_fields, vs.sample_data
    FROM variants v
    JOIN variant_samples vs ON v.variant_key = vs.variant_key
    WHERE v.count = 1
    ''')
    
    unique_variants = cursor.fetchall()
    
    # Счетчик уникальных вариантов для каждого образца
    sample_unique_counts = {}
    
    for variant in unique_variants:
        variant_key, chrom, pos, id, ref, alt, sample_name, qual, filter_fields, info_fields, format_fields, sample_data = variant
        
        # Увеличиваем счетчик уникальных вариантов для образца
        if sample_name in sample_unique_counts:
            sample_unique_counts[sample_name] += 1
        else:
            sample_unique_counts[sample_name] = 1
        
        # Записываем вариант в файл
        variant_line = [
            chrom, pos, id, ref, alt, qual, filter_fields, info_fields, format_fields, sample_data
        ]
        
        with open(f"{path_output}/{sample_name}_unique.vcf", 'a') as f:
            f.write("\t".join(variant_line) + "\n")
    
    # Обновляем количество уникальных вариантов в базе данных
    for sample_name, count in sample_unique_counts.items():
        cursor.execute('''
        UPDATE samples SET count_unique = ? WHERE sample_name = ?
        ''', (count, sample_name))
    
    conn.commit()
    conn.close()

def write_low_frequency_variants( folder_output,data_base, max_match=2):
    """
    Записывает варианты с низкой частотой встречаемости в файл.
    """
    print("Записываем низкочастотные snp")
    check_folder(folder_output)
    conn = sqlite3.connect(f"{folder_output}/{data_base}")
    cursor = conn.cursor()
    
    # Получаем варианты с низкой частотой встречаемости
    cursor.execute('''
    SELECT v.variant_key, v.chrom, v.pos, v.ref, v.alt, v.count
    FROM variants v
    WHERE v.count <= ?
    ''', (max_match,))
    
    low_frequency_variants = cursor.fetchall()
    low_frequency_SNP = len(low_frequency_variants)
    
    with open(f"{folder_output}/low_frequency_variants.txt", 'w') as f:
        f.write("Chrom\tPos\tRef\tAlt\tSample_Quals\tSample_Depths\tSample_Ref_Alt_Depths\n")
        
        for variant in low_frequency_variants:
            variant_key, chrom, pos, ref, alt, count = variant
            
            # Получаем данные о качестве, глубине и т.д. для этого варианта
            cursor.execute('''
            SELECT sample_name, qual, depth, ref_alt_depth
            FROM variant_samples
            WHERE variant_key = ?
            ''', (variant_key,))
            
            variant_samples = cursor.fetchall()
            
            sample_quals = []
            sample_depths = []
            sample_ref_alt_depths = []
            
            for vs in variant_samples:
                _, qual, depth, ref_alt_depth = vs
                sample_quals.append(qual)
                sample_depths.append(depth if depth else "")
                sample_ref_alt_depths.append(ref_alt_depth if ref_alt_depth else "")
            
            sample_quals_str = "/".join(sample_quals)
            sample_depths_str = "/".join(sample_depths)
            sample_ref_alt_depths_str = "/".join(sample_ref_alt_depths)
            
            variant_line = f"{chrom}\t{pos}\t{ref}\t{alt}\t{sample_quals_str}\t{sample_depths_str}\t{sample_ref_alt_depths_str}\n"
            f.write(variant_line)
    conn.close()
    return low_frequency_SNP


def apply_filter(cursor,filters_applied,filter_name, filter_query, filter_value):
    """
    Применяет фильтр и сохраняет данные о количестве отфильтрованных вариантов.
    """
    # Получаем количество вариантов до фильтрации
    cursor.execute('''
    SELECT sample_name, COUNT(*)
    FROM filtered_variants
    GROUP BY sample_name
    ''')
    counts_before = {row[0]: row[1] for row in cursor.fetchall()}
    
    # Применяем фильтр
    cursor.execute(filter_query, (filter_value,))
    filters_applied.append(filter_name)
    
    # Получаем количество вариантов после фильтрации
    cursor.execute('''
    SELECT sample_name, COUNT(*)
    FROM filtered_variants
    GROUP BY sample_name
    ''')
    counts_after = {row[0]: row[1] for row in cursor.fetchall()}
    
    # Сохраняем данные о фильтре
    for sample_name in counts_before:
        count_before = counts_before.get(sample_name, 0)
        count_after = counts_after.get(sample_name, 0)
        filtered_count = count_before - count_after  # Количество отфильтрованных вариантов
        
        cursor.execute('''
        INSERT INTO filters (sample_name, filter_name, count)
        VALUES (?, ?, ?)
        ''', (sample_name, filter_name, filtered_count))
    
def filter_variants(folder_output, data_base, min_depth=0, max_depth=999999, min_qual=0, min_alt_depth=0, min_vaf=0, max_vaf=1):
    """
    Фильтрует варианты по заданным параметрам и сохраняет данные о фильтрах в таблицу filters.
    """
    print("Фильтруем данные по заданным параметрам")
    conn = sqlite3.connect(f"{folder_output}/{data_base}")
    cursor = conn.cursor()
    
    # Удаляем временную таблицу, если она существует
    cursor.execute('''
    DROP TABLE IF EXISTS filtered_variants
    ''')
    
    # Создаем временную таблицу для отфильтрованных вариантов
    cursor.execute('''
    CREATE TABLE filtered_variants AS
    SELECT v.variant_key, v.chrom, v.pos, v.ref, v.alt, v.count,
           vs.sample_name, vs.qual, vs.filter_fields, vs.info_fields, 
           vs.format_fields, vs.sample_data, vs.depth, vs.ref_alt_depth, vs.alt_depth, vs.VAF
    FROM variants v
    JOIN variant_samples vs ON v.variant_key = vs.variant_key
    WHERE v.count = 1
    ''')
    
    # Применяем фильтры и сохраняем данные о фильтрах
    filters_applied = []
    
    # Применяем фильтры
    if min_depth > 0:
        apply_filter(
            cursor,filters_applied,
            filter_name=f"depth>={min_depth}",
            filter_query='''
            DELETE FROM filtered_variants
            WHERE depth IS NULL OR CAST(depth AS INTEGER) < ?
            ''',
            filter_value=min_depth
        )
    
    if max_depth < 999999:
        apply_filter(
            cursor,filters_applied,
            filter_name=f"depth<={max_depth}",
            filter_query='''
            DELETE FROM filtered_variants
            WHERE depth IS NULL OR CAST(depth AS INTEGER) > ?
            ''',
            filter_value=max_depth
        )
    
    if min_qual > 0:
        apply_filter(
            cursor,filters_applied,
            filter_name=f"qual>={min_qual}",
            filter_query='''
            DELETE FROM filtered_variants
            WHERE qual IS NULL OR CAST(qual AS FLOAT) < ?
            ''',
            filter_value=min_qual
        )
    
    if min_alt_depth > 0:
        apply_filter(
            cursor,filters_applied,
            filter_name=f"alt_depth>={min_alt_depth}",
            filter_query='''
            DELETE FROM filtered_variants
            WHERE alt_depth IS NULL OR CAST(alt_depth AS FLOAT) < ?
            ''',
            filter_value=min_alt_depth
        )
    
    if min_vaf > 0:
        apply_filter(
            cursor,filters_applied,
            filter_name=f"min_vaf>={min_vaf}",
            filter_query='''
            DELETE FROM filtered_variants
            WHERE VAF IS NULL OR CAST(VAF AS FLOAT) < ?
            ''',
            filter_value=min_vaf
        )
    
    if max_vaf < 1:
        apply_filter(
            cursor,filters_applied,
            filter_name=f"max_vaf<={max_vaf}",
            filter_query='''
            DELETE FROM filtered_variants
            WHERE VAF IS NULL OR CAST(VAF AS FLOAT) > ?
            ''',
            filter_value=max_vaf
        )
    
    conn.commit()
    conn.close()
    
    return filters_applied


def filter_variants_min_vaf(folder_output, data_base, min_vaf=0):
    """
    Фильтрует варианты по min_vaf и сохраняет данные.
    """
    print("Фильтруем данные по min_vaf")
    conn = sqlite3.connect(f"{folder_output}/{data_base}")
    cursor = conn.cursor()
    
    
    # Создаем временную таблицу для отфильтрованных вариантов
    cursor.execute('''
    CREATE TABLE filtered_variants_min_vaf AS
    SELECT *
    FROM filtered_variants
    ''')
    
    # Применяем фильтры и сохраняем данные о фильтрах
    filters_applied = []

    
    if min_vaf > 0:
        cursor.execute('''
        DELETE FROM filtered_variants_min_vaf
        WHERE VAF IS NULL OR CAST(VAF AS FLOAT) < ?
        ''', (min_vaf,))
        filters_applied.append(f"min_vaf>={min_vaf}")
    
    conn.commit()
    conn.close()
    
    return filters_applied


def filter_variants_max_vaf(folder_output, data_base, max_vaf=99999):
    """
    Фильтрует варианты по max_vaf и сохраняет данные.
    """
    print("Фильтруем данные по max_vaf")
    conn = sqlite3.connect(f"{folder_output}/{data_base}")
    cursor = conn.cursor()
    
    
    # Создаем временную таблицу для отфильтрованных вариантов
    cursor.execute('''
    CREATE TABLE filtered_variants_max_vaf AS
    SELECT *
    FROM filtered_variants
    ''')
    
    # Применяем фильтры и сохраняем данные о фильтрах
    filters_applied = []

    
    if max_vaf < 1:
        cursor.execute('''
        DELETE FROM filtered_variants_max_vaf
        WHERE VAF IS NULL OR CAST(VAF AS FLOAT) > ?
        ''', (max_vaf,))
        filters_applied.append(f"max_vaf<={max_vaf}")
        # Сохраняем данные о фильтре
    
    conn.commit()
    conn.close()
    
    return filters_applied



def write_unique_variants_with_filter(name_database, folder_output,filter_params ,table_name, add_coment_in_output_file):
    """
    Записывает отфильтрованные уникальные варианты в файлы.
    """
    print(f"Записываем данные с пройденым фильтром {filter_params} для кадого образца")
    conn = sqlite3.connect(f"{folder_output}/{name_database}")
    cursor = conn.cursor()
    
    # Создаем папку для отфильтрованных вариантов
    path_output = f"{folder_output}/unique_variants_{filter_params}"
    check_folder(path_output)
    
    # Получаем все образцы
    cursor.execute('''
    SELECT DISTINCT sample_name
    FROM filtered_variants
    ''')
    
    samples = cursor.fetchall()
    
    sample_data = ["filter"]
    sample_count_unique = [f"{filter_params}"]
    
    for sample in samples:
        sample_name = sample[0]
        
        # Создаем файл для образца
        with open(f"{path_output}/{sample_name}_unique_{filter_params}.txt", 'w') as f:
            if add_coment_in_output_file:
                # Добавляем заголовки из исходного VCF файла
                cursor.execute('''
                SELECT header_line
                FROM vcf_headers
                WHERE sample_name = ?
                ''', (sample_name,))
                
                headers = cursor.fetchall()
                for header in headers[:-1]:  # Все заголовки кроме последнего
                    f.write(header[0])
            
            # Добавляем заголовок таблицы
            f.write("Chrom\tPos\tRef\tAlt\tQual\tDepth\tRef_Alt_Depth\n")
        
        sample_data.append(sample_name)
        
        # Получаем варианты для образца
        cursor.execute(f'''
        SELECT chrom, pos, ref, alt, qual, depth, ref_alt_depth
        FROM {table_name}
        WHERE sample_name = ?
        ''', (sample_name,))
        
        variants = cursor.fetchall()
        sample_count_unique.append(len(variants))
        
        # Записываем варианты
        with open(f"{path_output}/{sample_name}_unique_{filter_params}.txt", 'a') as f:
            for variant in variants:
                chrom, pos, ref, alt, qual, depth, ref_alt_depth = variant
                variant_line = [chrom, pos, ref, alt, qual, depth, ref_alt_depth]
                f.write("\t".join(map(str, variant_line)) + "\n")
    
    # Записываем статистику
    with open(f"{folder_output}/Statistics_variant_with_parametr.txt", "a", encoding="utf-8") as file:
        file.write("\t".join(sample_data) + "\n")
        file.write("\t".join(map(str, sample_count_unique)) + "\n")
    conn.commit()
    conn.close()


def table_statistica(folder_output, name_database, LOW_FRAQUENCY_SNP, max_match_variant, filters_applied):
    """
    Создает таблицу статистики по вариантам.
    """
    print("Создаем файл с данными статистики")
    conn = sqlite3.connect(f"{folder_output}/{name_database}")
    cursor = conn.cursor()
    
    # Получаем все образцы
    cursor.execute('''
    SELECT sample_name, count_unique
    FROM samples
    ''')
    
    samples = cursor.fetchall()
    
    # Получаем все фильтры из базы данных
    cursor.execute('''
    SELECT DISTINCT filter_name
    FROM filters
    ''')
    
    filters = cursor.fetchall()
    filter_names = [filter[0] for filter in filters]
    
    # Упорядочиваем фильтры в соответствии с filters_applied
    ordered_filters = []
    for filter_applied in filters_applied:
        if filter_applied in filter_names:
            ordered_filters.append(filter_applied)
    
    # Создаем заголовок таблицы
    table = [["Sample", "unique_values(before_filter)", "unique_values(after_filter)"]]
    table[0].extend(ordered_filters)  # Добавляем упорядоченные названия фильтров в заголовок
    
    # Заполняем таблицу данными
    for sample in samples:
        sample_name, count_unique = sample
        row = [sample_name, count_unique]
        
        # Получаем количество вариантов после применения всех фильтров
        cursor.execute('''
        SELECT COUNT(*)
        FROM filtered_variants
        WHERE sample_name = ?
        ''', (sample_name,))
        
        filtered_count = cursor.fetchone()[0]
        row.append(filtered_count)  # Уникальные значения после фильтрации
        
        # Добавляем данные по каждому фильтру в порядке ordered_filters
        for filter_name in ordered_filters:
            cursor.execute('''
            SELECT count
            FROM filters
            WHERE sample_name = ? AND filter_name = ?
            ''', (sample_name, filter_name))
            
            result = cursor.fetchone()
            filter_count = result[0] if result else 0
            row.append(filter_count)  # Количество отфильтрованных вариантов для данного фильтра
        
        table.append(row)
    
    # Записываем таблицу в файл
    with open(f"{folder_output}/Statistics_variant.txt", "w", encoding="utf-8") as file:
        for row in table:
            line = "\t".join(map(str, row))
            file.write(line + "\n")
        
        # Добавляем дополнительную информацию
        additional_line = f"Number of detected SNP occurring in no more than {max_match_variant} samples: {LOW_FRAQUENCY_SNP} (before filtration)"
        file.write(additional_line + "\n")
    
    conn.close()


def visualization_variant_frequencies(folder_output, name_database):
    """
    Визуализирует частоты встречаемости вариантов.
    """
    print("Создаем картинки рапределения частот snp")
    conn = sqlite3.connect(f"{folder_output}/{name_database}")
    
    cursor = conn.cursor()
    
    # Получаем общее количество образцов
    cursor.execute('''
    SELECT COUNT(*) FROM samples
    ''')
    total_samples = cursor.fetchone()[0]
    
    # Получаем частоты встречаемости вариантов
    cursor.execute('''
    SELECT count FROM variants
    ''')
    
    counts = cursor.fetchall()
    frequencies = [count[0] / total_samples for count in counts]
    
    # Гистограмма
    plt.figure(figsize=(10, 6))
    plt.hist(frequencies, bins=10, edgecolor='black')
    plt.title("Распределение частот встречаемости вариантов")
    plt.xlabel("Частота встречаемости")
    plt.ylabel("Количество вариантов")
    plt.grid(True)
    plt.savefig(f"{folder_output}/histogram_variant_frequencies.png", dpi=300, bbox_inches='tight')
    
    # Density plot
    plt.figure(figsize=(10, 6))
    sns.kdeplot(frequencies, fill=True)
    plt.title("Плотность распределения частот встречаемости вариантов")
    plt.xlabel("Частота встречаемости")
    plt.ylabel("Плотность")
    plt.grid(True)
    plt.savefig(f"{folder_output}/density_plot_variant_frequencies.png", dpi=300, bbox_inches='tight')
    conn.close()


def drop_all_tables(folder_output,db_name):
    """
    Удаляет все таблицы в базе данных SQLite.
    """
    conn = sqlite3.connect(f"{folder_output}/{db_name}")
    cursor = conn.cursor()

    # Получаем список всех таблиц в базе данных
    cursor.execute("SELECT name FROM sqlite_master WHERE type='table'")
    tables = cursor.fetchall()

    # Удаляем каждую таблицу
    for table in tables:
        table_name = table[0]
        cursor.execute(f"DROP TABLE IF EXISTS {table_name}")
        print(f"Таблица {table_name} удалена.")

    conn.commit()
    conn.close()


if __name__ == "__main__":
    vcf_files_path = "/mnt/8Tb_new/arkom/test_task-2/"  # Путь к директории с VCF файлами
    folder_output = "/mnt/8Tb_new/arkom/results_task-3"  # Папка для сохранения результатов
    name_database = "vcf_data.db"
    max_match_variant = 2       # Порог для низкой частоты встречаемости SNP
    add_coment_in_output_file = True # Добавлять ли комментарии из исходных VCF файлов в выходные файлы


    vcf_files = get_vcf_files(vcf_files_path)

    setup_database(folder_output,name_database)

    parse_vcf_files(vcf_files, vcf_files_path, folder_output,name_database)

    write_unique_variants( folder_output, add_coment_in_output_file,name_database)

    LOW_FRAQUENCY_SNP = write_low_frequency_variants( folder_output,name_database, max_match_variant)

    filters_applied = filter_variants(
                folder_output,
                name_database,
                min_depth=10,
                min_qual=30,
                min_alt_depth=5,
                min_vaf=0.2
            )

    table_statistica(folder_output,name_database, LOW_FRAQUENCY_SNP, max_match_variant,filters_applied)


    filters_applied_1 = filter_variants_min_vaf(
                folder_output,
                name_database,
                min_vaf=0.65
            )

    filters_applied_2 = filter_variants_max_vaf(
                folder_output,
                name_database,
                max_vaf=0.3
            )

    write_unique_variants_with_filter(name_database, folder_output,"max_vaf=0.3" ,"filtered_variants_max_vaf", add_coment_in_output_file)
    write_unique_variants_with_filter(name_database, folder_output,"min_vaf=0.65" ,"filtered_variants_min_vaf", add_coment_in_output_file)


    visualization_variant_frequencies(folder_output,name_database)

    drop_all_tables(folder_output,name_database)

