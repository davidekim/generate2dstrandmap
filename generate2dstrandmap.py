#!/usr/bin/env python3
import os, sys
import numpy as np
import re
import copy
import argparse
#from cairosvg import svg2png

import pyrosetta
from pyrosetta import rosetta

fsizetermini = 46  # font size for N and C termini
fsize = 24   # font size for resnums in circles
hblen = 30   # length of hbond dashed lines
hbwidth = 5  # width of hbond dashed lines
hbondxshift = 14 # amount of left or right shift of hbond for don or acc
arrowlen = 60 # backbone strand direction arrow length
linewidth = 8 # backbone line width
termlen = 10 # extra len to backbone line at N and C ends
shearstrandopacity = "0.4" # transparency of the shear strand if one exists
spacer = 20 # spacer between circles
radius = 30 # radius of circles
rbulge = 20  # radius of bulge circles
circlestrokewidth = 3 # width of circle line
mainlinecolor = "deepskyblue" # backbone line color
bulgecolor = "red"  # color of beta bulges
h3_10color = "yellow"  # color of 3-10 helices
color_G = 'lime'  # color of gly in strands (leave blank to skip)

# init
resnums_3_10 = [] # 3-10 helix resnums
# if there are not enough hbonds to determine orientation, you can manually set them
# an example is the domain of 3fp9A which only has 1 hbond between parallel strands 3 and 5
# so the fix would be 3,5,'p' set by using 3-5-p as an arg  
fixed_orientations = [] # [[strandnum_i,strandnum_j,'p or a']] 

find_3_10 = False
hbs = {} # hbonds NH - O
hbsacc = {} # hbonds O - NH (acc key)
aastr = "" # AA sequence
ssstr = "" # secondary structure
abegostr = "" # abego
shearstrandindex = -1
strands_top_to_bottom_order = []
strands = [] # ordered from N->C
pleat = {} # key: xcoord val: pleat
skip_aa = []
add_bulges = []
add_Es = []
rm_Es = []

def init():
    global hbs,hbsacc,aastr,ssstr,abegostr,shearstrandindex,strands_top_to_bottom_order,strands,pleat
    hbs = {} # hbonds NH - O
    hbsacc = {} # hbonds O - NH (acc key)
    aastr = "" # AA sequence
    ssstr = "" # secondary structure
    abegostr = "" # abego
    shearstrandindex = -1
    strands_top_to_bottom_order = []
    strands = [] # ordered from N->C
    pleat = {} # key: xcoord val: pleat

def score_pose(p):
    scorefxn_tal_name="beta_nov15"
    scorefxn = rosetta.core.scoring.ScoreFunctionFactory.create_score_function(scorefxn_tal_name)
    energymethodoptions=scorefxn.energy_method_options()
    energymethodoptions.hbond_options().decompose_bb_hb_into_pair_energies(True)
    scorefxn.set_energy_method_options( energymethodoptions );
    return scorefxn(p)

# get backbone hbonds
def find_hbonds(p,ssstr):
    global hbs,hbsacc,resnums_3_10
    # initialize hbs and hbsacc
    for i,aa in enumerate(ssstr):
        iplus = i+1
        if iplus not in hbs:
            hbs[iplus] = []
        if iplus not in hbsacc:
            hbsacc[iplus] = []
    x=score_pose(p)
    p.update_residue_neighbors();
    hbond_set = rosetta.core.scoring.hbonds.HBondSet()
    rosetta.core.scoring.hbonds.fill_hbond_set(p, False, hbond_set)
    hbond_set.setup_for_residue_pair_energies(p , False, False);
    for i_hbond in range(hbond_set.nhbonds() ):
        test_hbond=hbond_set.hbond(i_hbond+1)
        if (test_hbond.acc_atm_is_backbone() and test_hbond.don_hatm_is_backbone()):
            accRes=test_hbond.acc_res()
            donRes=test_hbond.don_res()
            atomAname=p.residue(test_hbond.acc_res()).atom_name(test_hbond.acc_atm())
            atomBname=p.residue(test_hbond.don_res()).atom_name(test_hbond.don_hatm())
            if ssstr[int(donRes)-1] == 'E' and  ssstr[int(accRes)-1] == 'E':
                hbsacc[int(accRes)].append(int(donRes))
                hbs[int(donRes)].append(int(accRes))
            elif ssstr[int(donRes)-1] == 'E' and int(accRes) in resnums_3_10:
                hbsacc[int(accRes)].append(int(donRes))
                hbs[int(donRes)].append(int(accRes))
            elif ssstr[int(accRes)-1] == 'E' and int(donRes) in resnums_3_10:
                hbsacc[int(accRes)].append(int(donRes))
                hbs[int(donRes)].append(int(accRes))

def save_pleat():
    xcoords = []
    for s in strands:
        for d in s['dat']:
            if not d['bp']['abego'] == 'A' and d['resnum'] not in add_bulges:
                xcoords.append(int(d['x']))
    xcoords = list(dict.fromkeys(xcoords))
    xcoords.sort()
    for i,x in enumerate(xcoords):
        if i%2 == 0:
            pleat[x] = 1
        else:
            pleat[x] = 0

def pleat_color(x):
    if x in pleat and pleat[x]:
        return 'lavender'
    else:
        return 'white'
    
def xy(resnum,s):
    for d in s['dat']:
        if d['resnum'] == resnum:
            return [d['x'],s['y']]
    return []
        
def resnums(s):
    n = []
    for d in s['dat']:
        n.append(d['resnum'])
    return n

def direction(s):
    if s['dat'][0]['resnum'] < s['dat'][1]['resnum']:
        return 'r'
    else:
        return 'l'

def in_fixed_orientations(si, sj):
    for o in fixed_orientations:
        xin = si['n']
        xjn = sj['n']
        if si['n'] in o and sj['n'] in o:
            return o[2]
    return ''

def orientation(si,sj):
    for i in resnums(si):
        for j in resnums(sj):
            fixed = in_fixed_orientations(si,sj)
            if fixed:
                return fixed
            if j in hbs[i] and i in hbs[j]:
                return 'a'
            if i in hbs[j]:
                for x in hbs[i]:
                    if x+2 == j:
                        return 'p'
        for i in resnums(sj):
            for j in resnums(si):
                if j in hbs[i] and i in hbs[j]:
                    return 'a'
                if i in hbs[j]:
                    for x in hbs[i]:
                        if x+2 == j:
                            return 'p'

    return ''

def initialize_strands():
    # determine direction/ordering/pairing of strands
    # and initialize x and y coords
    sinit = [0]
    i = 0
    strands[0]['y'] = ((radius*2)+hblen)*(len(strands)+2)
    pairs = {}
    while len(sinit) < len(strands):
        thisi = i
        for j, sj in enumerate(strands):
            if j not in sinit:
                o = orientation(strands[i], strands[j])
                idir = direction(strands[i])
                jdir = direction(strands[j])
                if o == 'a':
                    if idir == 'r':
                        if jdir == 'r':
                            strands[j]['dat'].reverse()
                    elif jdir == 'l':
                        strands[j]['dat'].reverse()
                elif o == 'p':
                    if idir == 'r':
                        if jdir == 'l':
                            strands[j]['dat'].reverse()
                    elif jdir == 'r':
                        strands[j]['dat'].reverse()
                if o:
                    if strands[i]['bottomstrand'] >= 0:
                        strands[i]['topstrand'] = j
                        strands[j]['bottomstrand'] = i
                    else:
                        strands[i]['bottomstrand'] = j
                        strands[j]['topstrand'] = i
                    # initialize y
                    strands[j]['y'] = strands[sinit[-1]]['y'] - (radius*2 + hblen)
                    # initialize x
                    xpos = 0
                    for k, pos in enumerate(strands[i]['dat']):
                        if strands[i]['dat'][k]['bp']['abego'] != 'A' and strands[i]['dat'][k]['resnum'] not in resnums_3_10 and strands[i]['dat'][k]['resnum'] not in add_bulges:
                            if k-1 >= 0 and (strands[i]['dat'][k-1]['bp']['abego'] == 'A' or strands[i]['dat'][k-1]['resnum'] in resnums_3_10) or strands[i]['dat'][k-1]['resnum'] in add_bulges:
                                strands[i]['dat'][k]['x'] = xpos + radius + spacer/2 #bulge
                            else:
                                strands[i]['dat'][k]['x'] = xpos + radius*2 + spacer
                        else: # beta bulge
                            strands[i]['dat'][k]['x'] = xpos + radius + spacer/2 #bulge
                        xpos = strands[i]['dat'][k]['x']
                    xpos = 0
                    for k, pos in enumerate(strands[j]['dat']):
                        if strands[j]['dat'][k]['bp']['abego'] != 'A' and strands[j]['dat'][k]['resnum'] not in resnums_3_10 and strands[j]['dat'][k]['resnum'] not in add_bulges:
                            strands[j]['dat'][k]['x'] = xpos + radius*2 + spacer
                            if k-1 >= 0 and (strands[j]['dat'][k-1]['bp']['abego'] == 'A' or strands[j]['dat'][k-1]['resnum'] in resnums_3_10) or strands[j]['dat'][k-1]['resnum'] in add_bulges:
                                strands[j]['dat'][k]['x'] = xpos + radius + spacer/2 #bulge
                        else: # beta bulge
                            strands[j]['dat'][k]['x'] = xpos + radius + spacer/2 #bulge
                        xpos = strands[j]['dat'][k]['x']
                    sinit.append(j)
                    pairs[f'{i} {j}'] = 1
                    i = j
                    break
        if thisi == i:
            i = i + 1
    # add (shear) strand to bottom value of the last strand if it's paired with another strand (to complete a barrel for example)
    for j, sj in enumerate(strands):
        if sinit[-1] == j or f'{i} {j}' in pairs or f'{j} {i}' in pairs:
            continue
        if orientation(strands[sinit[-1]], strands[j]):
            strands[sinit[-1]]['bottomstrand'] = j
            break
    global strands_top_to_bottom_order
    strands_top_to_bottom_order = sinit
                        

def svg_backbone():
    svgstr = ""
    # draw horizontal lines
    strandresnums = []
    halfarrowlen = int(arrowlen/2)
    overhang = radius+termlen
    for i,strand in enumerate(strands):
        for d in strand['dat']:
            strandresnums.append(d['resnum'])
        y = strand['y']
        x1 = strand['dat'][0]['x']
        x2 = strand['dat'][-1]['x']
        opacity = strand['opacity']
        svgstr = svgstr + f'<line x1="{x1-overhang}" y1="{y}" x2="{x2+overhang}" y2="{y}" style="stroke:{mainlinecolor};stroke-width:{linewidth};opacity:{opacity}" />' + "\n"
        # strand N->C direction arrows
        n = strand['n']
        if direction(strand) == 'r':
            cx = x2+overhang
            offset = 8
            if n > 9:
                offset = 0
            svgstr = svgstr + f'<polygon points="{cx},{y+halfarrowlen},{cx},{y-halfarrowlen},{cx+arrowlen},{y}" style="fill:{mainlinecolor};stroke:{mainlinecolor};stroke-width:1;opacity:{opacity}" />' + "\n"  
            svgstr = svgstr + f'<text font-size="{fsize}" x="{cx+(arrowlen/2)-fsize+offset}" y="{y+(fsize/3)}" opacity="{opacity}" font-weight="bold" fill="black">{n}</text>' + "\n"
        else:
            offset = 6
            if n > 9:
                offset = 10
            cx = x1-overhang
            svgstr = svgstr + f'<polygon points="{cx},{y+halfarrowlen},{cx},{y-halfarrowlen},{cx-arrowlen},{y}" style="fill:{mainlinecolor};stroke:{mainlinecolor};stroke-width:1;opacity:{opacity}" />' + "\n"  
            svgstr = svgstr + f'<text font-size="{fsize}" x="{cx-(arrowlen/2)+(fsize/2)-offset}" y="{y+(fsize/3)}" opacity="{opacity}" font-weight="bold" fill="black">{n}</text>' + "\n"
    strandresnums = list(dict.fromkeys(strandresnums)) # remove duplicates
    strandresnums.sort()
    Nxy = xy(strandresnums[0],strands[0])
    if direction(strands[0]) == 'r':
        Nxy[0] = Nxy[0]-overhang-fsizetermini
    else:
        Nxy[0] = Nxy[0]+overhang+fsizetermini
    svgstr = svgstr + f'<text font-size="{fsizetermini}" x="{Nxy[0]}" y="{Nxy[1]+int(fsizetermini/3)}" font-weight="bold" fill="black">N</text>' + "\n"
    for s in strands:
        Cxy = xy(strandresnums[-1],s)
        if len(Cxy) == 2:
            if direction(s) == 'r':
                Cxy[0] = Cxy[0]+overhang+arrowlen+fsizetermini/3
            else:
                Cxy[0] = Cxy[0]-overhang-arrowlen-fsizetermini
            svgstr = svgstr + f'<text font-size="{fsizetermini}" x="{Cxy[0]}" y="{Cxy[1]+int(fsizetermini/3)}" font-weight="bold" fill="black">C</text>' + "\n"
            break
    if shearstrandindex >= 0:
        Nxy = xy(strandresnums[0],strands[shearstrandindex])
        if direction(strands[shearstrandindex]) == 'r':
            Nxy[0] = Nxy[0]-overhang-fsizetermini #arrowlen
        else:
            Nxy[0] = Nxy[0]+overhang+fsizetermini/3 #arrowlen
        opacity = strands[shearstrandindex]['opacity'] 
        svgstr = svgstr + f'<text font-size="{fsizetermini}" x="{Nxy[0]}" y="{Nxy[1]+int(fsizetermini/3)}" opacity="{opacity}" font-weight="bold" fill="black">N</text>' + "\n"
    return svgstr

def svg_circles():
    svgstr = ""
    for i,strand in enumerate(strands):
        opacity = strand['opacity']
        bulgestr = ""
        groupstr = f'<g id="strand{i+1}">'   
        for pos in strand['dat']:
            cy = strand['y']
            cr = radius
            resnum = pos['resnum']
            fill = pleat_color(pos['x'])
            small = False
            add_small_G_circle = False
            cx = pos['x']
            fx = pos['x']-fsize/2
            fy = cy+fsize/3
            if resnum < 10:
                fx = fx + fsize/3
            if resnum > 99:
                fx = fx - fsize/4
            if color_G and aastr[resnum-1] == 'G':
                fill = color_G
            if pos['bp']['abego'] == 'A' or pos['resnum'] in add_bulges:
                fill = bulgecolor
                cr = rbulge
                small = True
            if resnum in resnums_3_10:
                fill = h3_10color
                cr = rbulge
                small = True
                if color_G and aastr[resnum-1] == 'G':
                    add_small_G_circle = True
            if small:
                bulgestr = bulgestr + f'<circle cx="{cx}" cy="{cy}" r="{cr}" stroke="black" fill="{fill}" stroke-opacity="{opacity}" stroke-width="{circlestrokewidth}"/>' + "\n"  
                if add_small_G_circle:
                    bulgestr = bulgestr + f'<circle cx="{cx}" cy="{cy}" r="{cr-8}" stroke="black" fill="{color_G}" stroke-opacity="0" stroke-width="0"/>' + "\n"
                bulgestr = bulgestr + f'<text font-size="{fsize}" x="{fx}" y="{fy}" opacity="{opacity}" font-weight="bold" fill="black">{resnum}</text>' + "\n"                    
            else:
                groupstr = groupstr + f'<circle cx="{cx}" cy="{cy}" r="{cr}" stroke="black" fill="{fill}" stroke-opacity="{opacity}" stroke-width="{circlestrokewidth}"/>' + "\n"   
                groupstr = groupstr + f'<text font-size="{fsize}" x="{fx}" y="{fy}" opacity="{opacity}" font-weight="bold" fill="black">{resnum}</text>' + "\n"        
        groupstr = groupstr + bulgestr + "</g>"
        svgstr = svgstr + groupstr
    return svgstr

def hb_exists(hb,hbonds):
    hb.sort
    if hb in hbonds:
        return True
    return False

def svg_hbonds():
    svghbonds = ""
    for i,strand in enumerate(strands):
        top = strands[strand['topstrand']]
        for dat in strand['dat']:
            res = dat['resnum']
            resxy = [dat['x'],strand['y']]
            xshift = -hbondxshift
            if direction(strand) == 'l':
                xshift = hbondxshift
            if res in hbs:
                for acc in hbs[res]:
                    if acc in resnums(top):
                        accxy = xy(acc,top)
                        if abs(resxy[0] - accxy[0])<= (radius+rbulge)*2:
                            hxshift = xshift
                            if orientation(strand,top) == 'p':
                                hxshift = 0
                            svghbonds = svghbonds + f'<line x1="{resxy[0]+hxshift}" y1="{resxy[1]}" x2="{accxy[0]+hxshift}" y2="{accxy[1]}" style="stroke:rgb(255,0,0);stroke-width:{hbwidth};stroke-dasharray:{hbwidth}" />' + "\n"
            if res in hbsacc:
                for don in hbsacc[res]:
                    if don in resnums(top):
                        donxy = xy(don,top)
                        if abs(donxy[0] - resxy[0])<= (radius+rbulge)*2:
                            hxshift = -xshift
                            if orientation(strand,top) == 'p':
                                hxshift = 0  
                            svghbonds = svghbonds + f'<line x1="{donxy[0]+hxshift}" y1="{donxy[1]}" x2="{resxy[0]+hxshift}" y2="{resxy[1]}" style="stroke:rgb(255,0,0);stroke-width:{hbwidth};stroke-dasharray:{hbwidth}" />' + "\n"  
    return svghbonds
        
def add_shearstrand():
    global strands_top_to_bottom_order
    bottomstrand = strands_top_to_bottom_order[-1]
    if strands[strands_top_to_bottom_order[-1]]['bottomstrand'] >= 0:
        shearstrand = copy.deepcopy(strands[strands[strands_top_to_bottom_order[-1]]['bottomstrand']])
        shearstrand['bottomstrand'] = -1
        shearstrand['topstrand'] = strands_top_to_bottom_order[-1]
        strands[strands_top_to_bottom_order[-1]]['bottomstrand'] = strands_top_to_bottom_order[-1]+1
        shearstrand['opacity'] = shearstrandopacity
        shearstrand['y'] = strands[strands_top_to_bottom_order[-1]]['y']-(radius*2+hblen)
        o = orientation(strands[strands_top_to_bottom_order[-1]],shearstrand)
        bottomdir = direction(strands[strands_top_to_bottom_order[-1]])
        sheardir = direction(shearstrand)
        reversed = 0
        if o == 'p':
            if bottomdir == 'r' and sheardir == 'l':
                shearstrand['dat'].reverse()
                reversed = 1
            elif bottomdir == 'l' and sheardir == 'r':
                shearstrand['dat'].reverse()
                reversed = 1
        else:
            if bottomdir == 'r' and sheardir == 'r':
                shearstrand['dat'].reverse()
                reversed = 1
            elif bottomdir == 'l' and sheardir == 'l':
                shearstrand['dat'].reverse()
                reversed = 1
        if reversed:
            x = []
            for i,d in enumerate(shearstrand['dat']):
                x.append(d['x'])
            startx = x[-1]
            for i,d in enumerate(shearstrand['dat']):
                shearstrand['dat'][i]['x'] = startx+x[0]-x[i] 
        strands.append(shearstrand)
        global shearstrandindex
        shearstrandindex = len(strands)-1
        strands_top_to_bottom_order.append(shearstrandindex)

def add_shear_info():
    svgstr = ""
    spacetoline = hblen/2
    if shearstrandindex > 0:
        x1 = strands[strands_top_to_bottom_order[-1]]['dat'][0]['x']
        if direction(strands[strands_top_to_bottom_order[-1]]) == 'l':
            x1 = strands[strands_top_to_bottom_order[-1]]['dat'][-1]['x']
        y1 = strands[strands_top_to_bottom_order[-1]]['y']
        y2 = strands[strands_top_to_bottom_order[0]]['y']
        x2 = strands[strands_top_to_bottom_order[0]]['dat'][0]['x']
        x1gt = False
        if x1 > x2:
            x1gt = True
            x1 = strands[strands_top_to_bottom_order[-1]]['dat'][-1]['x']
            if direction(strands[strands_top_to_bottom_order[-1]]) == 'l':
                x1 = strands[strands_top_to_bottom_order[-1]]['dat'][0]['x']
            x2 = strands[strands_top_to_bottom_order[0]]['dat'][-1]['x']
        for s in strands_top_to_bottom_order[1:]:
            for d in strands[s]['dat']:
                if d['x'] == x2:
                    y2 = strands[s]['y']
        strandn = len(strands)-1

        y1 = y1 - (radius + spacetoline)
        y1_2 = y1 - (radius + hblen/2)
        shear = int(abs((x2-x1)/(radius*2+spacer)))
        has3_10 = False
        for s in strands: # check if there is a 3-10 helix within the shear span
            for d in s['dat']:
                if d['resnum'] in resnums_3_10:  
                    if x1gt and d['x'] <= x1 and d['x'] >= x2:
                        has3_10 = True
                        break
                    elif not x1gt and d['x'] >= x1 and d['x'] <= x2:
                        has3_10 = True
                        break
        if has3_10:
            shear = shear-1
        print(f'n: {strandn} S: {shear}')
        fsize = 30
        svgstr = svgstr + f'<line x1="{x1}" y1="{y1}" x2="{x1}" y2="{y1_2+1}" style="stroke:black;stroke-width:3;opacity:1" />' + "\n"
        svgstr = svgstr + f'<line x1="{x1}" y1="{y1_2}" x2="{x2}" y2="{y1_2}" style="stroke:black;stroke-width:3;opacity:1" />' + "\n"
        #svgstr = svgstr + f'<line x1="{x2}" y1="{y1_2}" x2="{x2}" y2="{y2+radius+spacetoline}" style="stroke:black;stroke-width:3;opacity:1" />' + "\n"
        svgstr = svgstr + f'<line x1="{x2}" y1="{y1_2+1}" x2="{x2}" y2="{y1}" style="stroke:black;stroke-width:3;opacity:1" />' + "\n"
        svgstr = svgstr + f'<text font-size="{fsize}" x="{(x1+(x2-x1)/2)-3*fsize}" y="{y1_2+fsize}" opacity="1" font-weight="bold" fill="black">S = {shear} (n = {strandn})</text>' + "\n"        
    return svgstr

def pair_strands():
    for i, si_index in enumerate(strands_top_to_bottom_order):
        if i == len(strands_top_to_bottom_order)-1:
            break
        sj_index = strands_top_to_bottom_order[i+1]
        si = strands[si_index]
        sj = strands[sj_index]
        o = orientation(si,sj)
        shift = None
        for k, ires in enumerate(resnums(si)):
            if abegostr[ires-1] == 'A' or ires in resnums_3_10 or ires in add_bulges: # skip bulge or 3-10
                continue
            for l, jres in enumerate(resnums(sj)):
                if abegostr[jres-1] == 'A' or jres in resnums_3_10 or jres in add_bulges : # skip bulge or 3-10
                    continue
                if jres in hbs[ires]:
                    if o == 'p':
                        shift = int(xy(jres,sj)[0]-xy(ires-1,si)[0])
                    else:
                        shift = int(xy(jres,sj)[0]-xy(ires,si)[0])
                elif ires in hbs[jres]:
                    if o == 'p':
                        shift = int(xy(jres,sj)[0]-xy(ires+1,si)[0])
                    else:
                        shift = int(xy(jres,sj)[0]-xy(ires,si)[0])
                if shift is not None:
                    break
            if shift is not None:
                break
        if shift is None:
            shift = 0
        for j in range(i,-1,-1) :
            x = strands_top_to_bottom_order[j]
            for k, dat in enumerate(strands[x]['dat']):
                strands[x]['dat'][k]['x'] = strands[x]['dat'][k]['x']+shift

def zero_x():
    global strands
    minx = 999999
    overhang = radius+termlen
    for s in strands:
        for d in s['dat']:
            if d['x'] < minx:
                minx = d['x']
    for i,s in enumerate(strands):
        for j,d in enumerate(s['dat']):
            if minx < 0:
                strands[i]['dat'][j]['x'] = strands[i]['dat'][j]['x']+abs(minx)+overhang+arrowlen+20+fsizetermini
            else:
                strands[i]['dat'][j]['x'] = strands[i]['dat'][j]['x']-minx+overhang+arrowlen+20+fsizetermini

def dimensions():
    global strands
    minx = 999999
    maxx = -999999
    miny = 999999
    maxy = -999999
    overhang = radius+termlen
    for s in strands:
        for d in s['dat']:
            if d['x'] < minx:
                minx = d['x']    
            if d['x'] > maxx:
                maxx = d['x']
            if s['y'] < miny:
                miny = s['y']    
            if s['y'] > maxy:
                maxy = s['y']
    return [maxx-minx+((fsizetermini+overhang+arrowlen)*2)+20,(maxy-miny)+radius*6]

def print_hdons():
    hdons = {}
    for don in hbs:
        for acc in hbs[don]:
            if don not in hdons:
                hdons[don] = []
            else:
                hdons[don].append(acc)
    for acc in hbsacc:
        for don in hbsacc[acc]:
            if don not in hdons:
                hdons[don] = []
            else:
                hdons[don].append(acc)
    for don in hdons:
        for acc in hdons[don]:
            print(f'don:{don} acc:{acc}')

def print_3_10():
    if resnums_3_10:
        resnums_3_10_str = [str(i) for i in resnums_3_10]
        r310str = ", ".join(resnums_3_10_str)
        print(f'3-10 helix: {r310str}')


parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('--strand_orientation', type=str, help='Manually set strand i to j orientation a or p. Example: 1-2-p,2-3-a')
parser.add_argument('--3_10', type=str, help='Manually set 3-10 helix by giving the starting resnum and length. Example: 52,3')
parser.add_argument('--find_3_10', type=bool, default=False, help='Try to identify a 3-10 helix.')
parser.add_argument('--add_hbonds', type=str, help='Manually set hbonds (don:acc) since they can be missed. Example: 60:37,94:98')
parser.add_argument('--add_E', type=str, help='Manually set secondary structure assignment to E. Example: 5,6,8-20,9')
parser.add_argument('--rm_E', type=str, help='Manually remove E secondary structure assignment. Example: 5,6,8-20,9')
parser.add_argument('--add_bulges', type=str, help='Manually set bulges since they can be missed. Example: 60,94')
parser.add_argument('--skip_aa', type=str, help='Skip strands with residues that are not part of the sheet. Example: 70,81-90,101')
parser.add_argument('pdbs', nargs=argparse.REMAINDER)
args = vars(parser.parse_args())
exit = False
if args['strand_orientation']:
    ofix = args['strand_orientation']
    for s in ofix.split(','):
        fix = s.split('-')
        if len(fix) == 3 and (fix[2] == 'p' or fix[2] == 'a'):
            fixed_orientations.append([int(fix[0]),int(fix[1]),fix[2]])
        else:
            exit = True
if args['3_10']:
    r310 = args['3_10'].split(',')
    if len(r310) != 2:
        exit = True
    else:
        for i in range(int(r310[0]),int(r310[0])+int(r310[1])):
            resnums_3_10.append(i)
if args['find_3_10']:
    find_3_10 = True
if args['skip_aa']:
    saa = args['skip_aa']
    for s in saa.split(','):
        skipaa = s.split('-')
        if len(skipaa) == 1:
            skip_aa.append(int(skipaa[0]))
        elif len(skipaa) == 2:
            for i in range(int(skipaa[0]),int(skipaa[1])+1):
                skip_aa.append(i)
        else:
            exit = True
if args['add_E']:
    addE = args['add_E']
    for s in addE.split(','):
        addEs = s.split('-')
        if len(addEs) == 1:
            add_Es.append(int(addEs[0]))
        elif len(addEs) == 2:
            for i in range(int(addEs[0]),int(addEs[1])+1):
                add_Es.append(i)
        else:
            exit = True
if args['rm_E']:
    rmE = args['rm_E']
    for s in rmE.split(','):
        rmEs = s.split('-')
        if len(rmEs) == 1:
            rm_Es.append(int(rmEs[0]))
        elif len(rmEs) == 2:
            for i in range(int(rmEs[0]),int(rmEs[1])+1):
                rm_Es.append(i)
        else:
            exit = True
if args['add_bulges']:
    bs = args['add_bulges']
    for b in bs.split(','):
        add_bulges.append(int(b))
if args['pdbs']:
    pdbs = args['pdbs']
else:
    exit = True
if exit:
    parser.print_help(sys.stderr)
    sys.exit(1)

import pyrosetta
from pyrosetta import rosetta

pyrosetta.init('-mute all -out::file::renumber_pdb -in::ignore_unrecognized_res -in::ignore_waters')

for i,pdb in enumerate(pdbs):
    init()

    if args['add_hbonds']:
        addhbs = args['add_hbonds'].split(',')
        if len(addhbs) > 0:
            for hbstr in addhbs:
                donacc = hbstr.split(':')
                if len(donacc) == 2:
                    if int(donacc[0]) not in hbs:
                        hbs[int(donacc[0])] = [int(donacc[1])]
                    else:
                        hbs[int(donacc[0])].append(int(donacc[1]))
                    if int(donacc[1]) not in hbsacc:
                        hbsacc[int(donacc[1])] = [int(donacc[0])]
                    else:
                        hbsacc[int(donacc[1])].append(int(donacc[0]))
                else:
                    parser.print_help(sys.stderr)
                    sys.exit(1)
    
    p = pyrosetta.pose_from_file(pdb)

    # get DSSP
    DSSP = pyrosetta.rosetta.core.scoring.dssp.Dssp(p)
    ssstr = DSSP.get_dssp_secstruct()
    aastr = p.sequence()
    # get ABEGO
    blankstr = ""
    abegostr = blankstr.join(pyrosetta.rosetta.core.sequence.ABEGOManager().get_symbols(p))

    tmpssstr = ssstr
    if add_Es:
        for res in add_Es:
            tmpssstr = tmpssstr[:res-1] + 'E' + tmpssstr[res:]
    if rm_Es:
        for res in rm_Es:
            tmpssstr = tmpssstr[:res-1] + 'L' + tmpssstr[res:]

    # get strand data
    if find_3_10:
        regex = re.compile(r"[E]+(H{3})?[E]+") # try to find a 3aa 3-10 helix?
    elif len(resnums_3_10) > 0: # fix 3-10 helix?
        for res in resnums_3_10:
            tmpssstr = tmpssstr[:res-1] + 'H' + tmpssstr[res:]
        len310 = len(resnums_3_10)
        regex = re.compile(r"[E]+(H{"+str(len310)+"})?[E]+") # 3-10 helix?
    else:
        regex = re.compile(r"[E]+")
    for match in re.finditer(regex, tmpssstr):
        if match.end()-match.start() < 3: # strands must have at least 3 E residues
            continue
        dat = []
        resnums3_10 = []
        if 'H' in match.group(0) and match.group(1):
            resnums3_10 = list(range(match.start(1)+1,match.end(1)+1))
        skipstrand = False
        for pos in range(match.start(),match.end()):
            pos = pos+1
            is3_10 = False
            if pos in resnums3_10:
                is3_10 = True
                if len(resnums_3_10) < 3:
                    resnums_3_10.append(pos)
            if pos in skip_aa:
                skipstrand = True
            dat.append({ 'resnum':pos, 'is3_10':is3_10, 'bp':{'aa':aastr[pos-1],'ss':ssstr[pos-1],'abego':abegostr[pos-1]}, 'x':0, 'r':0 })
        if not skipstrand:
            strands.append({ 'n':len(strands)+1,'dat':dat, 'y':0, 'topstrand':-1, 'bottomstrand':-1, 'opacity':1 })
    print(pdb)
    rulerstr = " "*len(aastr)
    for i,resn in enumerate(aastr,start=1):
        if i%10 == 0:
            rulerstr = rulerstr[:i-len(str(i))]+str(i)+rulerstr[i:]
    print(rulerstr)
    print(aastr)
    print(ssstr)
    if tmpssstr != ssstr:
        print(tmpssstr)
    print(abegostr)

    hbonds = find_hbonds(p,tmpssstr)
    initialize_strands()
    add_shearstrand()
    pair_strands()
    zero_x()
    save_pleat()
    print_hdons()
    print_3_10()
    svgstr = add_shear_info()
    svgstr = svgstr + svg_backbone() + svg_hbonds() + svg_circles()
    dims = dimensions()
    argsstring = " ".join(sys.argv)
    argsstring = argsstring.replace("--",'\\-\\-')
    svg = "<!--\n"+argsstring+"\n-->\n"
    svg = svg + f'<svg xmlns="http://www.w3.org/2000/svg" xmlns:xlink= "http://www.w3.org/1999/xlink" width="{dims[0]}" height="{dims[1]+radius*2}">'+"\n"
    svg = svg + f'<g transform="scale(1.0) rotate(0) translate(0,0)">'+svgstr+'</g></svg>'+"\n"

    f = open(pdb+'.2dstrandmap.svg','w')
    f.write(svg)
    f.close()

    #png = pdb+'.2dstrandmap.png'
    #svg2png(bytestring=svg,write_to=png)

    print()
