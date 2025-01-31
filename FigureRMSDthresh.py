import matplotlib.pyplot as plt


ref = ('8VTW', '#1/1A:_2375 #1/1A:_2376 #1/1A:_2377 #1/1A:_2378')


rmsd6_ref, rmsd6_all, rmsd4_ref = [], [], []

with open("GNRA_pairwise_RMSD.csv") as inp:
    for line in inp:
        rmsd = float(line.strip().split()[-1])
        if rmsd == 10.0:
            continue
        if ref[0] in line and ref[1] in line:
            rmsd6_ref.append(rmsd)
        rmsd6_all.append(rmsd)

with open("GNRA_pairwise_RMSD4.csv") as inp:
    for line in inp:
        rmsd = float(line.strip().split()[-1])
        if rmsd == 10.0 or rmsd == 0.0:
            continue
        rmsd4_ref.append(rmsd)

print(len(rmsd4_ref),len(rmsd6_all),len(rmsd6_ref))

plt.hist(rmsd6_all,bins=100)
plt.xlabel('RMSD', fontsize = 14)
plt.yticks([],[])
plt.xlim([-0.05,3.05])
plt.xticks([0.0,0.5,1.0,1.5,2.0,2.5,3.0],[0.0,0.5,1.0,1.5,2.0,2.5,3.0])
plt.savefig("RMSD_size6_all.svg")
plt.show()


plt.hist(rmsd6_ref,bins=20)
plt.xlabel('RMSD', fontsize = 14)
plt.yticks([],[])
plt.xlim([-0.05,3.05])
plt.xticks([0.0,0.5,1.0,1.5,2.0,2.5,3.0],[0.0,0.5,1.0,1.5,2.0,2.5,3.0])
plt.savefig("RMSD_size6_ref.svg")
plt.show()

plt.hist(rmsd4_ref,bins=20)
plt.xlabel('RMSD', fontsize = 14)
plt.yticks([],[])
plt.xlim([-0.05,3.05])
plt.xticks([0.0,0.5,1.0,1.5,2.0,2.5,3.0],[0.0,0.5,1.0,1.5,2.0,2.5,3.0])
plt.savefig("RMSD_size4_ref.svg")
plt.show()
