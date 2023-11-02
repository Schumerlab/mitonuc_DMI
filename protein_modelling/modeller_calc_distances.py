# Run in Pymol GUI:

cmd.fetch("6g2j")
cmd.load("input_files/xbir_NDUFS5.pdb")
cmd.align("xbir_NDUFS5","6g2j")
cmd.load("input_files/xbir_ND6.pdb")
cmd.align("xbir_ND6","6g2j")
cmd.load("input_files/Xbir_NDUFA13.pdb")


for idx in range(1,6):cmd.load("input_files/xbir_nd6-6g2j.B9999000%01d.pdb"%idx)
for idx in range(1,6):cmd.align("xbir_nd6-6g2j.B9999000%01d"%idx, "xbir_ND6")
for idx in range(1,6):cmd.color("cyan", "xbir_nd6-6g2j.B9999000%01d"%idx)
for idx in range(1,6):cmd.show("sticks", "xbir_nd6-6g2j.B9999000%01d and resi 7+45+46+52+119+120+122+126+135+147"%idx)

for idx in range(1,6):cmd.load("input_files/xbir_nd6-5lnk.B9999000%01d.pdb"%idx)
for idx in range(1,6):cmd.align("xbir_nd6-5lnk.B9999000%01d"%idx, "xbir_ND6")
for idx in range(1,6):cmd.color("magenta", "xbir_nd6-5lnk.B9999000%01d"%idx)
for idx in range(1,6):cmd.show("sticks", "xbir_nd6-5lnk.B9999000%01d and resi 7+45+46+52+119+120+122+126+135+147"%idx)

for idx in range(1,6):cmd.load("input_files/xbir_nd6-5ldw.B9999000%01d.pdb"%idx)
for idx in range(1,6):cmd.align("xbir_nd6-5ldw.B9999000%01d"%idx, "xbir_ND6")
for idx in range(1,6):cmd.color("yellow", "xbir_nd6-5ldw.B9999000%01d"%idx)
for idx in range(1,6):cmd.show("sticks", "xbir_nd6-5ldw.B9999000%01d and resi 7+45+46+52+119+120+122+126+135+147"%idx)

for idx in range(1,6):cmd.load("input_files/xbir_nd6-5xtc.B9999000%01d.pdb"%idx)
for idx in range(1,6):cmd.align("xbir_nd6-5xtc.B9999000%01d"%idx, "xbir_ND6")
for idx in range(1,6):cmd.color("green", "xbir_nd6-5xtc.B9999000%01d"%idx)
for idx in range(1,6):cmd.show("sticks", "xbir_nd6-5xtc.B9999000%01d and resi 7+45+46+52+119+120+122+126+135+147"%idx)


# open dist.txt for writing
# First, 31A in ndufs5
f=open('output_files/xbir_31M_dist.txt','w')

f.write(' '.join(['xbir_ND6','31M','119E',"%.3f\n"%cmd.distance('tmp','/xbir_NDUFS5///MET`31/CA','/xbir_ND6///GLU`119/CA')]))
f.write(' '.join(['xbir_ND6','31M','120A',"%.3f\n"%cmd.distance('tmp','/xbir_NDUFS5///MET`31/CA','/xbir_ND6///ALA`120/CA')]))
f.write(' '.join(['xbir_ND6','31M','122D',"%.3f\n"%cmd.distance('tmp','/xbir_NDUFS5///MET`31/CA','/xbir_ND6///ASP`122/CA')]))
f.write(' '.join(['xbir_ND6','31M','126I',"%.3f\n"%cmd.distance('tmp','/xbir_NDUFS5///MET`31/CA','/xbir_ND6///ILE`126/CA')]))
f.write(' '.join(['xbir_ND6','31M','135V',"%.3f\n"%cmd.distance('tmp','/xbir_NDUFS5///MET`31/CA','/xbir_ND6///VAL`135/CA')]))

for obj in cmd.get_object_list('xbir_nd6-*'):f.write(' '.join([obj,'31M','119E',"%.3f\n"%cmd.distance('tmp','/xbir_NDUFS5///MET`31/CA','/'+obj+'///GLU`119/CA')]))
for obj in cmd.get_object_list('xbir_nd6-*'):f.write(' '.join([obj,'31M','120A',"%.3f\n"%cmd.distance('tmp','/xbir_NDUFS5///MET`31/CA','/'+obj+'///ALA`120/CA')]))
for obj in cmd.get_object_list('xbir_nd6-*'):f.write(' '.join([obj,'31M','122D',"%.3f\n"%cmd.distance('tmp','/xbir_NDUFS5///MET`31/CA','/'+obj+'///ASP`122/CA')]))
for obj in cmd.get_object_list('xbir_nd6-*'):f.write(' '.join([obj,'31M','126I',"%.3f\n"%cmd.distance('tmp','/xbir_NDUFS5///MET`31/CA','/'+obj+'///ILE`126/CA')]))
for obj in cmd.get_object_list('xbir_nd6-*'):f.write(' '.join([obj,'31M','135V',"%.3f\n"%cmd.distance('tmp','/xbir_NDUFS5///MET`31/CA','/'+obj+'///VAL`135/CA')]))
# close the output file.
f.close()


# Then, 79Y in ndufa13 (note that distance is measured from tyrosine hydroxyl, not alpha carbon)
f=open('output_files/xbir_79Y_dist.txt','w')

f.write(' '.join(['xbir_ND6','79Y','119E',"%.3f\n"%cmd.distance('tmp','/Xbir_NDUFA13///TYR`79/OH','/xbir_ND6///GLU`119/CA')]))
f.write(' '.join(['xbir_ND6','79Y','120A',"%.3f\n"%cmd.distance('tmp','/Xbir_NDUFA13///TYR`79/OH','/xbir_ND6///ALA`120/CA')]))
f.write(' '.join(['xbir_ND6','79Y','122D',"%.3f\n"%cmd.distance('tmp','/Xbir_NDUFA13///TYR`79/OH','/xbir_ND6///ASP`122/CA')]))
f.write(' '.join(['xbir_ND6','79Y','126I',"%.3f\n"%cmd.distance('tmp','/Xbir_NDUFA13///TYR`79/OH','/xbir_ND6///ILE`126/CA')]))
f.write(' '.join(['xbir_ND6','79Y','135V',"%.3f\n"%cmd.distance('tmp','/Xbir_NDUFA13///TYR`79/OH','/xbir_ND6///VAL`135/CA')]))

for obj in cmd.get_object_list('xbir_nd6-*'):f.write(' '.join([obj,'79Y','119E',"%.3f\n"%cmd.distance('tmp','/Xbir_NDUFA13///TYR`79/OH','/'+obj+'///GLU`119/CA')]))
for obj in cmd.get_object_list('xbir_nd6-*'):f.write(' '.join([obj,'79Y','120A',"%.3f\n"%cmd.distance('tmp','/Xbir_NDUFA13///TYR`79/OH','/'+obj+'///ALA`120/CA')]))
for obj in cmd.get_object_list('xbir_nd6-*'):f.write(' '.join([obj,'79Y','122D',"%.3f\n"%cmd.distance('tmp','/Xbir_NDUFA13///TYR`79/OH','/'+obj+'///ASP`122/CA')]))
for obj in cmd.get_object_list('xbir_nd6-*'):f.write(' '.join([obj,'79Y','126I',"%.3f\n"%cmd.distance('tmp','/Xbir_NDUFA13///TYR`79/OH','/'+obj+'///ILE`126/CA')]))
for obj in cmd.get_object_list('xbir_nd6-*'):f.write(' '.join([obj,'79Y','135V',"%.3f\n"%cmd.distance('tmp','/Xbir_NDUFA13///TYR`79/OH','/'+obj+'///VAL`135/CA')]))
# close the output file.
f.close()

# Last, 140L in ndufa13 (measured from alpha carbon again)

f=open('output_files/xbir_140L_dist.txt','w')

f.write(' '.join(['xbir_ND6','140L','119E',"%.3f\n"%cmd.distance('tmp','/Xbir_NDUFA13///LEU`140/CA','/xbir_ND6///GLU`119/CA')]))
f.write(' '.join(['xbir_ND6','140L','120A',"%.3f\n"%cmd.distance('tmp','/Xbir_NDUFA13///LEU`140/CA','/xbir_ND6///ALA`120/CA')]))
f.write(' '.join(['xbir_ND6','140L','122D',"%.3f\n"%cmd.distance('tmp','/Xbir_NDUFA13///LEU`140/CA','/xbir_ND6///ASP`122/CA')]))
f.write(' '.join(['xbir_ND6','140L','126I',"%.3f\n"%cmd.distance('tmp','/Xbir_NDUFA13///LEU`140/CA','/xbir_ND6///ILE`126/CA')]))
f.write(' '.join(['xbir_ND6','140L','135V',"%.3f\n"%cmd.distance('tmp','/Xbir_NDUFA13///LEU`140/CA','/xbir_ND6///VAL`135/CA')]))

for obj in cmd.get_object_list('xbir_nd6-*'):f.write(' '.join([obj,'140L','119E',"%.3f\n"%cmd.distance('tmp','/Xbir_NDUFA13///LEU`140/CA','/'+obj+'///GLU`119/CA')]))
for obj in cmd.get_object_list('xbir_nd6-*'):f.write(' '.join([obj,'140L','120A',"%.3f\n"%cmd.distance('tmp','/Xbir_NDUFA13///LEU`140/CA','/'+obj+'///ALA`120/CA')]))
for obj in cmd.get_object_list('xbir_nd6-*'):f.write(' '.join([obj,'140L','122D',"%.3f\n"%cmd.distance('tmp','/Xbir_NDUFA13///LEU`140/CA','/'+obj+'///ASP`122/CA')]))
for obj in cmd.get_object_list('xbir_nd6-*'):f.write(' '.join([obj,'140L','126I',"%.3f\n"%cmd.distance('tmp','/Xbir_NDUFA13///LEU`140/CA','/'+obj+'///ILE`126/CA')]))
for obj in cmd.get_object_list('xbir_nd6-*'):f.write(' '.join([obj,'140L','135V',"%.3f\n"%cmd.distance('tmp','/Xbir_NDUFA13///LEU`140/CA','/'+obj+'///VAL`135/CA')]))
# close the output file.
f.close()
