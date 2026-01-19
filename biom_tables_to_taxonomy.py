
#If not already downoladed, you wil need the nodes, names and merged files from taxdump
# https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz
# cd /media/herve/DATA2/Scripts_and_comandlines
# usage python biom_tables_to_taxonomy.py

biom_file='/media/herve/10TB/Apogee/6_mock/15_analysis_bracken_v02/bracken_out/table.biom'  # Input file with taxID in first colum and named #otu 
taxonomy_file='/media/herve/10TB/Apogee/6_mock/15_analysis_bracken_v02/bracken_out/bracken_taxo.tsv'  # The output file
nodes_file='/media/herve/10TB/Apogee/1_data_base/Kraken/kraken_db/taxonomy/nodes.dmp'  
names_file='/media/herve/10TB/Apogee/1_data_base/Kraken/kraken_db/taxonomy/names.dmp'
merged_file='/media/herve/10TB/Apogee/1_data_base/Kraken/kraken_db/taxonomy/merged.dmp'

def expand_taxonomy_from_taxid(biom_file, taxonomy_file, nodes_file, names_file, merged_file):
        # Parse biom_file
        taxid_dict = dict()
        with open(biom_file, 'r') as f:
            for line in f:
                if line.startswith('#'):
                    continue
                line = line.rstrip()
                taxid = line.split('\t')[0]  # taxID in first field
                taxid_dict[taxid] = ''

        # Parse nodes.dmp file
        node_dict = dict()
        with open(nodes_file, 'r') as f:
            for line in f:
                line = line.rstrip()
                fields = line.split('\t')
                taxid = fields[0]
                parent = fields[2]
                rank = fields[4]

                node_dict[taxid] = (parent, rank)

        # Parse names.dmp file
        names_dict = dict()
        with open(names_file, 'r') as f:
            for line in f:
                if 'scientific name' in line:
                    fields = line.split('\t')
                    taxid = fields[0]
                    name = fields[2]
                    if taxid not in names_dict:
                        names_dict[taxid] = name
                    else:
                        continue

        # Parse meged.dmp
        with open(merged_file, 'r') as f:
            for line in f:
                fields = line.split('\t')
                old_taxid = fields[0]
                new_taxid = fields[2]
                if old_taxid in taxid_dict:
                    taxid_dict[new_taxid] = taxid_dict[old_taxid]
                    del taxid_dict[old_taxid]

        with open(taxonomy_file, 'w') as f:
            # Write header
            f.write('otu\tsuperkingdom\tphylum\tclade\torder\tfamily\tgenus\tspecies\n')

            # for taxid in taxid_list:
            for taxid, acc in taxid_dict.items():
                og_taxid = taxid
                taxo_list = list()
                while taxid != '1':
                    (parent, rank)=node_dict[taxid]
                    taxo_list.insert(0, (rank, names_dict[taxid]))
                    taxid = parent

                # Write taxonomy to file
                # Kingdom;Phylum;Class;Order;Family;Genus;Species
                # k__;p__;c__;o__;f__;g__;s__
                # unidentified
                taxo_dict = {'k': 'unidentified', 'p': 'unidentified', 'c': 'unidentified',
                             'o': 'unidentified', 'f': 'unidentified', 'g': 'unidentified', 's': 'unidentified'}

                for i in taxo_list:
                    if i[0] == 'species':
                        taxo_dict['s'] = i[1].replace(' ', '_')
                    elif i[0] == 'genus':
                        taxo_dict['g'] = i[1]
                    elif i[0] == 'family':
                        taxo_dict['f'] = i[1]
                    elif i[0] == 'order':
                        taxo_dict['o'] = i[1]
                    elif i[0] == 'clade':
                        taxo_dict['c'] = i[1]
                    elif i[0] == 'phylum':
                        taxo_dict['p'] = i[1]
                    elif i[0] == 'superkingdom':
                        taxo_dict['k'] = i[1]
                    else:
                        continue

                f.write('{}\t{}\t\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(
                    og_taxid, taxo_dict['k'], taxo_dict['p'], taxo_dict['c'],
                    taxo_dict['o'], taxo_dict['f'], taxo_dict['g'], taxo_dict['s']))


expand_taxonomy_from_taxid(biom_file, taxonomy_file, nodes_file, names_file, merged_file)
