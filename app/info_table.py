import pickle as pkl
import csv

with open('pigs.pkl', 'rb') as pklfile:
    pigs_clusters = pkl.load(pklfile)

with open('poultry.pkl', 'rb') as pklfile:
    poultry_clusters = pkl.load(pklfile)

with open('pigs.csv', 'wb') as csvfile:
    writer = csv.writer(csvfile)
    for key,value in pigs_clusters.iteritems():
        writer.writerow([key])
        writer.writerow([value['num_members']])
        for amr_name, amr_val in value['top_amrs']:
            writer.writerow([amr_name, amr_val/float(value['num_members'])*100])
        for country_name, country_val in value['predominant_country']:
            writer.writerow([country_name, country_val/float(value['num_members'])*100])
        # amr_perc = [(amr_name, amr_val/float(value['num_members'])*100) for amr_name, amr_val in value['top_amrs']]
        # country_perc = [(country_name, country_val/float(value['num_members'])*100) for country_name, country_val in value['predominant_country']]
        # writer.writerow([key,amr_perc,country_perc])


with open('poultry.csv', 'wb') as csvfile:
    writer = csv.writer(csvfile)
    for key,value in poultry_clusters.iteritems():
        writer.writerow([key])
        writer.writerow([value['num_members']])
        for amr_name, amr_val in value['top_amrs']:
            writer.writerow([amr_name, amr_val/float(value['num_members'])*100])
        for country_name, country_val in value['predominant_country']:
            writer.writerow([country_name, country_val/float(value['num_members'])*100])
