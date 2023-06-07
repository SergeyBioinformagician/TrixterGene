try:
    import random
    from Bio.Seq import MutableSeq
    from Bio import SeqIO
    import networkx as nx
    import pandas as pd
    import numpy as np
    import matplotlib.pyplot as plt
    import re
    from io import StringIO
except Exception as e:
    print('module connection error:', str(e))
    #exit()

# Объекты библиотеки TrixterGene
# Чтение таблиц


class Table:
    
    def table_creation(self,path,table_class):
        if table_class == "rep_annotation":
            self.df = pd.read_csv(path, sep= " ",header = None)
            self.df.columns = ['Smith_Waterman','%_substitutions',
                           '%_of_bases_opposite_a_gap','%_of_bases_opposite_a_gap_in_the_repeat_consensus',
                           'name_of_query_sequence','starting_position_of_match','ending_position_of_match',
                           'number_of_bases_in_query_sequence_past_the_ending_position_of_match',
                           'match_is_with_the_Complement_of_the_consensus_sequence_in_the_database',
                           'name_of_the_matching_interspersed_repeat',
                           'the_class_of_the_repeat',
                           'no._of_bases_in_(complement_of)_the_repeat_consensus_sequence...',
                           'starting_position_of_match_in_database_sequence',
                           'ending_position_of_match_in_database_sequence',
                           'ID','stars']
            self.stats = {col_name:val.describe() for col_name, val in self.df.iteritems()}
        elif table_class == "report":
            self.description = ''
            self.df_report = None
            desArr = []
            csvStr = ''
            with open(path,'r') as file:
                while True:
                    line = file.readline()
                    if line is None or line == '':
                        break
                    if line[0] == '#':
                        desArr.append(line)
                        continue
                    csvStr += line
            if len(desArr) > 0:
                csvStr = desArr.pop() + csvStr
            self.description = ''.join(desArr)
            csvStr = csvStr[1:]
            csvStr = csvStr.lstrip()
            self.df_report = pd.read_csv(StringIO(csvStr),sep= '\t')
            
            

class Rep_annotation(Table):
    
    
    def __init__(self, path_in, path_out = None):
        self.path_in = path_in
        self.path_out = path_out
        if path_out != None:
            self.to_dataframe(self.path_in,self.path_out)
            self.table_creation(self.path_out,"rep_annotation")
        else:
            self.table_creation(self.path_in,"rep_annotation")
        
              
    def parse_str(self,string):
        res = ''
        res = re.sub(" +", " ", string)
        return res
    
    
    def to_dataframe(self,path_in,path_out):
        
        try:    
            file1 = open(path_in, 'r')
        except IOError:
            print("Input file opening error")
      
        try:
            fileRes = open(path_out, 'w')
        except IOError:
            print("Output file opening error")
            
        file1.readline()
        file1.readline()
        file1.readline()
        
        while True:
            line = file1.readline()
            if not line:
                break
            if '*' not in line:
                line = line.rstrip() + ' ' + '\n'
            fileRes.write(self.parse_str(line).lstrip())
        fileRes.close()
        file1.close()


class Report(Table):
    
    def __init__(self,path):
        self.table_creation(path,'report')

# Чтение output результатов работы blast и работа с ними
class Blast_table:
    def __init__(self, file_path, my_columns=['qseqid', 'sseqid', 'pident', 'length', 'qlen', 'slen']):
        self.df = pd.read_csv(file_path, sep='\t', header=None, names=my_columns)
        
    def create_div_df(self):
        self.div_df = (self.df.query('qseqid != sseqid')).groupby('qseqid')['distance'].min().reset_index()
        self.div_df.rename(columns={'distance':'last_div_dist'}, inplace=True)

    def distance_col(self, model = None, my_model = False, **kwargs): # Реализует расчёт значений в столбце
    # с расстояниями в зависимости от модели 
        if my_model:
            kwargs['self_Blast_table'] = self
            result_series = self.df.apply(my_model, axis=1,**kwargs)
            self.df['distance'] = result_series

        pident = self.df['pident'].values
        length = self.df['length'].values
        qlen = self.df['qlen'].values
        slen = self.df['slen'].values

        if model == 'JC':
            D_JC = np.where((1 - ((pident*length)/np.maximum(qlen,slen))/100) < 0.75,
                            (-3/4)*np.log(1 + (-4/3)*(1 - ((pident*length)/np.maximum(qlen,slen))/100)),np.nan)

            self.df['distance'] = D_JC
        else:
            raise ValueError('Unknown model')
            
    
    def blast_dbscan(self,epsilon,mu = 1):
        def all_possible_clusters(no_reverse_repeats_df): # Функция превращает все возможные кластеры во множества
            return no_reverse_repeats_df.groupby('qseqid').apply(lambda x: set(x['qseqid']).union(set(x['sseqid']))).tolist()
        def sets_to_clusters(all_possible_clusters,connected_column): # Функция объединяет все мн-ва, имеющие пересечения
            all_possible_clusters = np.array(all_possible_clusters) # Входящий лист -> np array
            for name in connected_column: # Итерация всех вершин, имеющих ребро
                logical = np.array([name in group for group in all_possible_clusters]) # Проверка наличия элемента 
                # во всех мн-вах массива
                new_cluster = set().union(*all_possible_clusters[logical]) # Все мн-ва с общим элементов объед.
                all_possible_clusters = np.append(all_possible_clusters[~logical],new_cluster) # Срез массива, без слитых мн-в
                # + полученное объединение
            return all_possible_clusters # Теперь среди возможных кластеров нет пересечений

        def merge_sets_with_mu(set_array, all_names, mu): # Объединяет компоненты связности для mu != 1
            set_list = [set(x) for x in set_array]
            while True:
                merged = False
                for i in range(len(set_list)):
                    for j in range(i+1, len(set_list)):
                        common = set_list[i] & set_list[j]
                        if len(common) >= mu:
                            merged = True
                            set_list[i] |= set_list[j]
                            del set_list[j]
                            break
                    if merged:
                        break
                if not merged:
                    break
            return np.array([set(x) for x in set_list])



        nodes = self.df.qseqid.unique() # Имена всех локусов
        clustering_df = self.df.query('qseqid != sseqid') # Удаляем из таблицы сравнение локусов с самими собой
        if mu == 1:
            distances = clustering_df.distance < epsilon # Логический вектор для образования связей между парами локусов
            distances = distances.values
            actual_clusters_df = clustering_df[distances] # В кластеры попадут только те элементы, у которых есть хотя бы 1 пара
            connected = {name for name in actual_clusters_df.qseqid} # Мн-во связанных рёбрами вершин
            anomalies = {name for name in clustering_df[~distances].qseqid} # Мн-во пар без рёбер 
            anomalies -= connected # Из мн-ва пар без рёбер вычитаем связанные вершины, получаем аномалии
            connected_column = actual_clusters_df.qseqid.unique() # Массив связанных локусов
            names = np.sort(np.array([actual_clusters_df.qseqid, actual_clusters_df.sseqid]), axis=0) # Сортировка всех имён
            # в столбцах, чтобы можно было удалить случаи A   B; B   A методом drop_duplicates
            no_reverse_repeats_df = pd.DataFrame(names.T, columns=['qseqid', 'sseqid']).drop_duplicates() # Больше нет строк
            # повторяющих одну и ту же пару
            try:
                all_possible_c = all_possible_clusters(no_reverse_repeats_df)
                clusters = sets_to_clusters(all_possible_c,connected_column)
                clusters_dict = {"Cluster " + str(index + 1) : value for index, value in enumerate(clusters)}
                clusters_dict['Anomalies'] = anomalies
                self.clusters = clusters_dict,epsilon,mu

            except AttributeError: # Позволяет избежать ошибки в методе sets_to_clusters, когда остались лишь аномалии
                clusters_dict = {'Anomalies' : anomalies}
                self.clusters = clusters_dict,epsilon,mu
        else:
            distances = clustering_df.distance < epsilon # Логический вектор для образования связей между парами локусов
            distances = distances.values
            actual_clusters_df = clustering_df[distances] # В кластеры попадут только те элементы, у которых есть хотя бы 1 пара
            connected = {name for name in actual_clusters_df.qseqid} # Мн-во связанных рёбрами вершин
            anomalies = {name for name in clustering_df[~distances].qseqid} # Мн-во пар без рёбер 
            anomalies -= connected # Из мн-ва пар без рёбер вычитаем связанные вершины, получаем аномалии
            connected_column = actual_clusters_df.qseqid.unique() # Массив связанных локусов
            names = np.sort(np.array([actual_clusters_df.qseqid, actual_clusters_df.sseqid]), axis=0) # Сортировка всех имён
            # в столбцах, чтобы можно было удалить случаи A   B; B   A методом drop_duplicates
            no_reverse_repeats_df = pd.DataFrame(names.T, columns=['qseqid', 'sseqid']).drop_duplicates() # Больше нет строк
            # повторяющих одну и ту же пару
            try:
                all_possible_c = merge_sets_with_mu(no_reverse_repeats_df)
                clusters = sets_to_clusters(all_possible_c,connected_column,mu)
                clusters_dict = {"Cluster " + str(index + 1) : value for index, value in enumerate(clusters)}
                clusters_dict['Anomalies'] = anomalies
                self.clusters = clusters_dict,epsilon,mu

            except AttributeError: # Позволяет избежать ошибки в методе merge_sets_with_mu, когда остались лишь аномалии
                clusters_dict = {'Anomalies' : anomalies}
                self.clusters = clusters_dict,epsilon,mu
    
    def create_real_graph(self,cluster_name,gc_content_fasta = False):
        cluster_df = self.df[self.df['qseqid'].isin(self.clusters[0][cluster_name]) & self.df['sseqid'].isin(self.clusters[0][cluster_name])]
        filtered_df = cluster_df[cluster_df['distance'] <= self.clusters[1]]
        filtered_df = filtered_df[['qseqid','sseqid','distance']]
        filtered_df = filtered_df.query('qseqid != sseqid')
        if gc_content_fasta is not False:
            filtered_df = self.add_gc_content(gc_content_fasta,filtered_df)
        my_real_graph = Real_graph(filtered_df)
        return my_real_graph
    
    def add_gc_content(self,fasta_file_path, df):
        # Считываем FASTA файл в словарь
        sequences_dict = SeqIO.to_dict(SeqIO.parse(fasta_file_path, "fasta"))
        # Функция для вычисления GC-контента
        def gc_content(seq):
            return np.round((seq.count("G") + seq.count("C") + seq.count("c") + seq.count("g")) / len(seq), 2)
        # Применяем функцию gc_content к каждой последовательности в словаре
        gc_contents = {k: gc_content(str(v.seq)) for k, v in sequences_dict.items()}
        # Добавляем столбец gc_content1 и gc_content2 к датафрейму
        df['gc_content1'] = df['qseqid'].map(gc_contents)
        df['gc_content2'] = df['sseqid'].map(gc_contents)
        return df


# Парсинг таблицы

def get_all_classes(repeats_df):
    classes = repeats_df.the_class_of_the_repeat.unique()
    return classes
    
def get_all_chromosomes(repeats_df):
    chromosomes = repeats_df.name_of_query_sequence.unique()
    return chromosomes
    

# График частот длин повторов определённой длины

def plot_repeat_length(df, repeat_class, repeat_name=None, bin_size=100):
    if repeat_name is not None:
        df = df[df['name_of_the_matching_interspersed_repeat'] == repeat_name]
        repeat_class = df.iloc[0]['the_class_of_the_repeat'] # узнаем класс повтора по имени
        title = f"Length distribution of {repeat_name} repeat"
    else:
        df = df[df['the_class_of_the_repeat'] == repeat_class]
        title = f"Length distribution of {repeat_class} repeats"
    
    repeat_lengths = df['ending_position_of_match'] - df['starting_position_of_match']
    bins = range(0, int(repeat_lengths.max()), bin_size)
    plt.hist(repeat_lengths, bins=bins)
    plt.title(title)
    plt.xlabel("Repeat length (bp)")
    plt.ylabel("Frequency")
    plt.show()


# Объект для хранения данных о повторах определённого класса (классов) на заданных хромосомах

class Loci: # rep_classes - array, chromosomes - array
    def __init__(self, df, rep_classes, chromosomes, some_loci = False): 
        if not some_loci:
            self.rep_classes = rep_classes
            self.chromosomes = chromosomes
            self.initial_df = df
            self.filtered_df = self.initial_df[self.initial_df['the_class_of_the_repeat'].isin(rep_classes)]
            self.intervals = {chrom: self.filtered_df.loc[(self.filtered_df.name_of_query_sequence == chrom)] for chrom in self.chromosomes}
            # {имя хромосомы: таблица повторов на хромосоме}
            


def fasta_maker(loci_object,input_fasta,output_fasta):
    with open(output_fasta,"w") as my_fasta:
        count = 0
        with open(input_fasta) as handle:
            for record in SeqIO.parse(handle, "fasta"):
                for chrom in loci_object.intervals.keys():
                    if chrom in record.id:
                        for index, row in loci_object.intervals[chrom].iterrows():
                            my_fasta.write(f">{chrom},{row['the_class_of_the_repeat']},{row['name_of_the_matching_interspersed_repeat']},{row['starting_position_of_match']}:{row['ending_position_of_match']}\n")
                            my_fasta.write(str(record.seq[row.starting_position_of_match:row.ending_position_of_match])+'\n')
                        print(chrom)
                        
  

# Цикл для создания фаста файлов по элементам на хромосомах
def save_fasta_with_elem(path,chromosomes,elements,ref_gen_path):  
    for chromosome in chromosomes:
        for element in elements:
            current_locus = Loci(annotation.df,[element],[chromosome])
            path = path+str(chromosome)+'_'+element.replace('/', '-')+'.fasta'
            fasta_maker(current_locus,ref_gen_path,path)


# Исправление имён в fasta файлах
def fasta_name_gaps(input_fasta, output_fasta):
    with open(input_fasta, 'r') as infile, open(output_fasta, 'w') as outfile:
        for line in infile:
            if line.startswith('>'):
                new_name = line.rstrip().replace(' ', '_')
                outfile.write(new_name + '\n')
            else:
                outfile.write(line)
                

def blast_dbscan(df,epsilon,mu = 1):
    def all_possible_clusters(no_reverse_repeats_df): # Функция превращает все возможные кластеры во множества
        return no_reverse_repeats_df.groupby('qseqid').apply(lambda x: set(x['qseqid']).union(set(x['sseqid']))).tolist()
    def sets_to_clusters(all_possible_clusters,connected_column): # Функция объединяет все мн-ва, имеющие пересечения
        all_possible_clusters = np.array(all_possible_clusters) # Входящий лист -> np array
        for name in connected_column: # Итерация всех вершин, имеющих ребро
            logical = np.array([name in group for group in all_possible_clusters]) # Проверка наличия элемента 
            # во всех мн-вах массива
            new_cluster = set().union(*all_possible_clusters[logical]) # Все мн-ва с общим элементов объед.
            all_possible_clusters = np.append(all_possible_clusters[~logical],new_cluster) # Срез массива, без слитых мн-в
            # + полученное объединение
        return all_possible_clusters # Теперь среди возможных кластеров нет пересечений
    
    def merge_sets_with_mu(set_array, all_names, mu): # Объединяет компоненты связности для mu != 1
        set_list = [set(x) for x in set_array]
        while True:
            merged = False
            for i in range(len(set_list)):
                for j in range(i+1, len(set_list)):
                    common = set_list[i] & set_list[j]
                    if len(common) >= mu:
                        merged = True
                        set_list[i] |= set_list[j]
                        del set_list[j]
                        break
                if merged:
                    break
            if not merged:
                break
        return np.array([set(x) for x in set_list])

    

    nodes = df.qseqid.unique() # Имена всех локусов
    clustering_df = df.query('qseqid != sseqid') # Удаляем из таблицы сравнение локусов с самими собой
    if mu == 1:
        distances = clustering_df.JK < epsilon # Логический вектор для образования связей между парами локусов
        distances = distances.values
        actual_clusters_df = clustering_df[distances] # В кластеры попадут только те элементы, у которых есть хотя бы 1 пара
        connected = {name for name in actual_clusters_df.qseqid} # Мн-во связанных рёбрами вершин
        anomalies = {name for name in clustering_df[~distances].qseqid} # Мн-во пар без рёбер 
        anomalies -= connected # Из мн-ва пар без рёбер вычитаем связанные вершины, получаем аномалии
        connected_column = actual_clusters_df.qseqid.unique() # Массив связанных локусов
        names = np.sort(np.array([actual_clusters_df.qseqid, actual_clusters_df.sseqid]), axis=0) # Сортировка всех имён
        # в столбцах, чтобы можно было удалить случаи A   B; B   A методом drop_duplicates
        no_reverse_repeats_df = pd.DataFrame(names.T, columns=['qseqid', 'sseqid']).drop_duplicates() # Больше нет строк
        # повторяющих одну и ту же пару
        try:
            all_possible_c = all_possible_clusters(no_reverse_repeats_df)
            clusters = sets_to_clusters(all_possible_c,connected_column)
            clusters_dict = {"Cluster " + str(index + 1) : value for index, value in enumerate(clusters)}
            clusters_dict['Anomalies'] = anomalies
            return clusters_dict
        
        except AttributeError: # Позволяет избежать ошибки в методе sets_to_clusters, когда остались лишь аномалии
            clusters_dict = {'Anomalies' : anomalies}
            return clusters_dict
    else:
        distances = clustering_df.JK < epsilon # Логический вектор для образования связей между парами локусов
        distances = distances.values
        actual_clusters_df = clustering_df[distances] # В кластеры попадут только те элементы, у которых есть хотя бы 1 пара
        connected = {name for name in actual_clusters_df.qseqid} # Мн-во связанных рёбрами вершин
        anomalies = {name for name in clustering_df[~distances].qseqid} # Мн-во пар без рёбер 
        anomalies -= connected # Из мн-ва пар без рёбер вычитаем связанные вершины, получаем аномалии
        connected_column = actual_clusters_df.qseqid.unique() # Массив связанных локусов
        names = np.sort(np.array([actual_clusters_df.qseqid, actual_clusters_df.sseqid]), axis=0) # Сортировка всех имён
        # в столбцах, чтобы можно было удалить случаи A   B; B   A методом drop_duplicates
        no_reverse_repeats_df = pd.DataFrame(names.T, columns=['qseqid', 'sseqid']).drop_duplicates() # Больше нет строк
        # повторяющих одну и ту же пару
        try:
            all_possible_c = merge_sets_with_mu(no_reverse_repeats_df)
            clusters = sets_to_clusters(all_possible_c,connected_column,mu)
            clusters_dict = {"Cluster " + str(index + 1) : value for index, value in enumerate(clusters)}
            clusters_dict['Anomalies'] = anomalies
            return clusters_dict
        
        except AttributeError: # Позволяет избежать ошибки в методе merge_sets_with_mu, когда остались лишь аномалии
            clusters_dict = {'Anomalies' : anomalies}
            return clusters_dict
            
        
# Объект для хранения графового представления кластера с реальными мобильными элементами

class Real_graph:
    def __init__(self, mobile_graph):
        self.mobile_graph = mobile_graph
    def graph_creation(self):
        G = nx.Graph()
        for i, row in self.mobile_graph.iterrows():
            G.add_edge(row['qseqid'], row['sseqid'], weight=row['distance'])
        self.mobile_graph_network = G
    def draw_graph(self,with_labels = False):
        if 'mobile_graph_network' not in self.__dict__:
            print('you need to execute graph_creation()')
            return 
        G = self.mobile_graph_network
        pos = nx.spring_layout(G)
        if with_labels is False:
            # рисуем граф
            nx.draw(G, pos, with_labels=False, node_color='skyblue', node_size=10, font_size=1)    
        else:
            # определяем веса ребер в графе
            edge_labels = {(u, v): d['weight'] for u, v, d in G.edges(data=True)}
            # рисуем граф
            nx.draw(G, pos, with_labels=True, node_color='skyblue', node_size=10, font_size=1)
            nx.draw_networkx_edge_labels(G, pos, edge_labels=edge_labels, font_size=1)
        plt.show()
        

# Класс - модель мобильного элемента

class Mobile_element:
    def __init__(self,sequence,death_prob = 0.1,transp_prob = 0.9,my_number = None,parent_name = None,initial_name = 'hydra'):
        self.sequence = sequence.copy() # sequence = {участок:'его последовательность'}
        self.my_number = my_number # Какой я по счёту дочерний элемент
        self.parent_name = parent_name # имя родителя
        self.name = None
        self.initial_name = initial_name
        self.generate_name()
        self.number_of_children = 0 # Число детей у мобильного элемента
        self.age = 0 # мой возраст
        self.is_alive = True # жив ли я
        self.death_prob = death_prob # Ожидаемая продолжительность жизни от рождения
        self.transp_prob = transp_prob
        self.birthday = 0
        
    def mutate(self,my_law, *args): # Дай пользователю право задать свой эвол.закон
        new_sequence = my_law(self.sequence,*args)
        self.sequence = new_sequence
    
    def get_gc_content(self):
        res_str = ''
        for key in self.sequence:
            res_str += self.sequence[key]
        sum_G = res_str.count('G') + res_str.count('g')
        sum_C = res_str.count('C') + res_str.count('c')
        gc_content = (sum_G + sum_C) / len(res_str)
        return gc_content
        
    def aging(self):
        # Увеличение возраста элемента
        self.age += 1
        
    def death(self):
        # Функция убивающая мобильный элемент
        if random.random() <= self.death_prob:
            self.is_alive = False
            return 1 # для мониторинга числа погибших транспозонов
        else:
            return 0 # для мониторинга числа погибших транспозонов
        
    def die(self):
        self.is_alive = False
        return 1
            
    def reproduce(self,transp_prob,death_prob):
        if self.is_alive == True:
            dice = random.random()
            if dice <= self.transp_prob: # NB! У ребёнка может быть уже другая вероятность, которая
                # идёт входным параметром!
                tmp_var = Mobile_element(sequence = self.sequence,death_prob = death_prob,
                                         transp_prob = transp_prob,
                                        my_number = self.number_of_children+1,
                                         parent_name = self.name)
                self.number_of_children += 1
                return tmp_var
            elif dice > transp_prob:
                return None
        else:
            return None
        
        
    def generate_name(self):
        if self.parent_name is not None:
            parent_name = self.parent_name
            self.name = parent_name + ',' + str(self.my_number)
        else:
            self.name = self.initial_name #initial_name # Имя первого элемента
        
    def Jukes_Cantor(self,mutation_rate):
        for part in self.sequence.keys():
            mutable_seq = MutableSeq(self.sequence[part])
            for i, nucleotide in enumerate(self.sequence[part]):
                if random.random() < mutation_rate:
                    new_nucleotide = random.choice('ACTG'.replace(nucleotide, ''))
                    mutable_seq[i] = new_nucleotide
            self.sequence[part] = str(mutable_seq)
           
           
# Симулятор размножения


class Time:
    def __init__(self,initial_name,sequence,death_prob = 0.1,transp_prob = 0.8,mutation_rate = 0.05):
        self.hydra = Mobile_element(sequence,death_prob = death_prob,transp_prob = transp_prob,initial_name = initial_name)
        self.tik_tak = 0
        self.mobile_dict = dict()
        self.mobile_dict[self.hydra.name] = self.hydra
        self.mobile_graph_network = None
        self.death_prob = death_prob
        self.transp_prob = transp_prob
        self.mutation_rate = mutation_rate
        self.df_features = pd.DataFrame(pd.np.empty((0, 4)))
        self.df_features.columns = ['name','age','is_alive','gc_content'] # добавь gc состав
        self.mobile_graph = pd.DataFrame(pd.np.empty((0, 7))) 
        self.mobile_graph.columns = ['name1','name2','distance','is_alive1','is_alive2','gc_content1','gc_content2']

    def get_distance(self,name1,name2,weights=None, my_model = False, **kwargs):
        if my_model is True:
            distance = my_model(self.mobile_dict[name1],self.mobile_dict[name2],**kwargs)
            return distance
        arr_S = []
        arr_len = []
        for key in self.mobile_dict[name1].sequence.keys():
            np_arr1 = np.array(list(self.mobile_dict[name1].sequence[key]))
            np_arr2 = np.array(list(self.mobile_dict[name2].sequence[key]))
            np_diff = np_arr1 == np_arr2
            S = np_diff.sum() / len(np_diff)
            arr_S.append(S)
            arr_len.append(len(np_diff))
        np_arr_S = np.array(arr_S)
        arr_D = 1 - np_arr_S
        D_JC = np.where(arr_D < 0.75,(-3/4)*np.log(1 + (-4/3)*arr_D),1) 
        if weights is not None:
            pass
        else:
            weights = [length / sum(arr_len) for length in arr_len]
        np_weights = np.array(weights)
        D_JC_weights_product = np.dot(np_weights,D_JC)
        distance = D_JC_weights_product / sum(weights)
        return distance
        
    def expansion(self,transp_limit,time_limit = None,transposition_dynamics = None,death_dynamics = None,model_for_evolution = None):
        num_created = 1
        num_dead = 0
        while len(self.mobile_dict) < transp_limit:
            self.tik_tak += 1
            tmp_arr = []
            for elem in self.mobile_dict.keys():
                new_head = None
                if transposition_dynamics is None:
                    new_head = self.mobile_dict[elem].reproduce(self.transp_prob,self.death_prob)
                elif transposition_dynamics is not None:
                    new_head = self.mobile_dict[elem].reproduce(transposition_dynamics(self,self.mobile_dict[elem]),
                                                                self.death_prob)
                if new_head is not None:
                    num_created += 1
                    new_head.birthday = self.tik_tak
                    tmp_arr.append(new_head)
                self.mobile_dict[elem].aging()
                if model_for_evolution is None:
                    self.mobile_dict[elem].Jukes_Cantor(self.mutation_rate)
                else:
                    self.mobile_dict[elem].sequence = model_for_evolution(self,self.mobile_dict[elem])
                if death_dynamics is None:
                    num_dead += self.mobile_dict[elem].death()
                else:
                    death_chance = death_dynamics(self,self.mobile_dict[elem])
                    if random.random() <= death_chance:
                        num_dead += self.mobile_dict[elem].die()
            for elem in tmp_arr:
                self.mobile_dict[elem.name] = elem
                if len(self.mobile_dict) == transp_limit:
                    break
                
                
            if time_limit is not None:
                if time_limit <= self.tik_tak:
                    break

            tmp_arr = []
            if num_created <= num_dead:
                    break
    
    def common_ancestor(self,name1,name2):
        name1 = name1.split(',')
        name2 = name2.split(',')
        ancestor = ''
        if len(name1) >= len(name2):
            for i in range(len(name2)):
                if name1[i] == name2[i]:
                    ancestor += name1[i]
                    ancestor += ','
                else:
                    break
            return ancestor[:-1]
        else:
            for i in range(len(name1)):
                if name1[i] == name2[i]:
                    ancestor += name1[i]
                    ancestor += ','
                else:
                    break
            return ancestor[:-1]
        
    #def div_time(self,mobile_dict,name1,name2,tik_tak):
        #print('tik_tak = ' + str(tik_tak))
        #common_ancestor_res = self.common_ancestor(name1,name2)
        #print('common_ancestor_res = ' +common_ancestor_res)
        #mobile_dict_res = mobile_dict[common_ancestor_res]
        #print('type mobile_dict_res = ' + str(type(mobile_dict_res)))
        #birthday_res = mobile_dict_res.birthday
        #print(mobile_dict_res.age)
        #print('type birthday_res = '+ str(type(birthday_res)))
        #res = tik_tak - birthday_res
        #print(res)
        #return tik_tak - mobile_dict[self.common_ancestor(name1,name2)].birthday
    
    #def distance(self,mobile_dict,name1,name2,tik_tak):
        #return self.div_time(mobile_dict,name1,name2,tik_tak)*2*self.mutation_rate
    
    def graph_creation(self,epsilon):
        G = nx.Graph()
        mobile_graph_tmp = self.mobile_graph[self.mobile_graph['distance'] <= epsilon]
        for i, row in mobile_graph_tmp.iterrows():
            if i % 1000000 == 0:    
                print(i)
            G.add_edge(row['name1'], row['name2'], weight=row['distance'])
            node_color = 'skyblue' if row['is_alive1'] else 'red'
            node_data = {'color': node_color}
            G.add_node(row['name1'], **node_data)
            node_color2 = 'skyblue' if row['is_alive2'] else 'red'
            node_data2 = {'color': node_color2}
            G.add_node(row['name2'], **node_data2)
        self.mobile_graph_network = G
    
    def draw_graph(self,node_color = False):
        pos = nx.spring_layout(self.mobile_graph_network)
        if node_color is False:
            # рисуем граф
            node_colors = [self.mobile_graph_network.nodes[n]['color'] for n in self.mobile_graph_network.nodes()]
            nx.draw(self.mobile_graph_network, pos, with_labels=False, node_color=node_colors, node_size=10, font_size=1)
        else:
            nx.draw(self.mobile_graph_network, pos, with_labels=False, node_color='skyblue', node_size=10, font_size=1)
        plt.show()
    
    
    def count_alive_mobiles(self):
        alive = 0
        dead = 0
        for mobile in self.mobile_dict.values():
            if mobile.is_alive:
                alive += 1
            else:
                dead += 1
        return alive, dead
    
    def add_row_to_df_features(self,name,age,is_alive):
        gc_content = self.mobile_dict[name].get_gc_content()
        new_row = {'name':name,'age':age,'is_alive':is_alive,'gc_content':gc_content}
        self.df_features.loc[len(self.df_features)] = new_row

    def mobile_dict_to_df_features(self):
        self.df_features = self.df_features[0:0] # удаление всех строк в df_features
        for key,value in self.mobile_dict.items():
            self.add_row_to_df_features(name = key,age = value.age,is_alive = value.is_alive)
        self.df_features['number_of_edges'] = self.df_features.groupby('name')['name'].transform('count')
    
    def add_row_to_mobile_graph(self,name1,name2,distance):
        is_alive1 = self.mobile_dict[name1].is_alive
        is_alive2 = self.mobile_dict[name2].is_alive
        gc_content1 = self.mobile_dict[name1].get_gc_content()
        gc_content2 = self.mobile_dict[name2].get_gc_content()
        new_row = {'name1':name1,'name2':name2,'distance':distance,'is_alive1':is_alive1,'is_alive2':is_alive2,
                  'gc_content1':gc_content1,'gc_content2':gc_content2}
        self.mobile_graph.loc[len(self.mobile_graph)] = new_row
    
    def mobile_dict_to_mobile_graph(self,epsilon):
        self.mobile_graph = self.mobile_graph[0:0] # удаление всех строк в mobile_graph
        dict_for_check = {}
        for key,value in self.mobile_dict.items():
            for key2,value2 in self.mobile_dict.items():
                if key != key2:
                    if str(key + key2) in dict_for_check:
                        continue
                    distance_res = self.get_distance(key,key2)#self.distance(self.mobile_dict,key,key2,self.tik_tak)
                    self.add_row_to_mobile_graph(key,key2,distance_res)
                    dict_for_check[key2+key] = 1
        self.mobile_graph = self.mobile_graph[self.mobile_graph['distance'] <= epsilon]
        self.mobile_graph['number_of_edges'] = self.mobile_graph.groupby('name1')['name1'].transform('count')
