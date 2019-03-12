def sphereclash(Spherelist):

    #Iterates through list of spheres and performs a pairwise comparison to determine if two given spheres have the potential to clash in space
    noclash = None
    whichclash = []
    for sphere in Spherelist:
        for sphere2 in Spherelist:
            if (sphere[1] != sphere2[1]) == True and noclash != False:
                spherecenter = (((sphere2[1][0] - sphere[1][0])**2 + (sphere2[1][1] - sphere[1][1])**2 + (sphere2[1][2] - sphere[1][2])**2)**0.5)
                rmax = sphere[0] + sphere2[0]
                noclash = spherecenter > rmax
                listobject = []
                listobject.append(sphere)
                listobject.append(sphere2)
                listobject.append(noclash)
                whichclash.append(listobject)
            if (sphere[1] != sphere2[1]) == True and noclash == False:
                spherecenter = (((sphere2[1][0] - sphere[1][0])**2 + (sphere2[1][1] - sphere[1][1])**2 + (sphere2[1][2] - sphere[1][2])**2)**0.5)
                rmax = sphere[0] + sphere2[0]
                currentclash = spherecenter > rmax
                listobject = []
                listobject.append(sphere)
                listobject.append(sphere2)
                listobject.append(currentclash)
                whichclash.append(listobject)

    return noclash, whichclash
