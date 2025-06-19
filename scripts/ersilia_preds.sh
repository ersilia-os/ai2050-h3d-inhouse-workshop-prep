ersilia -v serve eos3804
ersilia -v run -i ../data/libraries/antiinfective_smiles.csv -o ../data/libraries/antiinfective_eos3804.csv
ersilia close

ersilia -v serve eos42ez
ersilia -v run -i ../data/libraries/antiinfective_smiles.csv -o ../data/libraries/antiinfective_eos42ez.csv
ersilia close

ersilia -v serve eos2db3
ersilia -v run -i ../data/libraries/antiinfective_smiles.csv -o ../data/libraries/antiinfective_eos2db3.csv
ersilia close

ersilia -v serve eos7d58
ersilia -v run -i ../data/libraries/antiinfective_smiles.csv -o ../data/libraries/antiinfective_eos7d58.csv
ersilia close

ersilia -v serve eos18ie
ersilia -v run -i ../data/libraries/antiinfective_smiles.csv -o ../data/libraries/antiinfective_eos18ie.csv
ersilia close


ersilia -v serve eos3804
ersilia -v run -i ../data/libraries/generalistic_smiles.csv -o ../data/libraries/generalistic_eos3804.csv
ersilia close

ersilia -v serve eos42ez
ersilia -v run -i ../data/libraries/generalistic_smiles.csv -o ../data/libraries/generalistic_eos42ez.csv
ersilia close

ersilia serve eos2db3
ersilia -v run -i ../data/libraries/generalistic_smiles.csv -o ../data/libraries/generalistic_eos2db3.csv
ersilia close

ersilia serve eos7d58
ersilia -v run -i ../data/libraries/generalistic_smiles.csv -o ../data/libraries/generalistic_eos7d58.csv
ersilia close

ersilia serve eos18ie
ersilia -v run -i ../data/libraries/generalistic_smiles.csv -o ../data/libraries/generalistic_eos18ie.csv
ersilia close