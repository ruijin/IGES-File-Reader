#!/usr/bin/env python
import io
import numpy as np

from read_IGES import *
from spline import nurbs2bezier

def exportNURBSSurface(surface):
    file_str = io.StringIO()
    
    degu = surface.M1
    degv = surface.M2

    #print(surface)

    # first convert nurbs to bezier
    # TODO: assume the knot sequence has knot multiplicty degree + 1 at both end
    #print(len(surface.control_points), len(surface.W))
    cpts = np.array([np.array([surface.control_points[i][0],
                      surface.control_points[i][1],
                      surface.control_points[i][2],
                      surface.W[i]]) for i in range(len(surface.W))])

    cpts = cpts.reshape([ surface.K2 + 1, surface.K1 + 1, 4])

    # 1rd direction
    b_knot_u, b_cpts = nurbs2bezier(surface.T1[1:-1], [cpts[:,i] for i in range(surface.K1 + 1)], degu)
    cpts = np.array(b_cpts)
    b_knot_v, b_cpts = nurbs2bezier(surface.T2[1:-1], [cpts[:,i] for i in range(surface.K2 + 1)], degv)
    cpts = np.array(b_cpts)
    
    uvbbox = [(surface.U0, surface.V0), (surface.U1, surface.V1)]

    total_bicubic_piece = 0
    
    # TODO: assume the nurbs curve is indeed in bezier form
    for piece_v in range((cpts.shape[0] - 1)//degv):
        v_base = b_knot_v[piece_v * degv]
        v_end = b_knot_v[piece_v * degv + degv]
        for piece_u in range((cpts.shape[1] - 1)//degu):
            u_base = b_knot_u[piece_u * degu]
            u_end = b_knot_u[piece_u * degu + degu]

            bcpts = np.zeros([degv + 1, degu + 1, 6])
            
            for i in range(degv + 1):
                for j in range(degu + 1):
                    idx_v = i + degv * piece_v
                    idx_u = j + degu * piece_u

                    p = cpts[idx_v][idx_u]

                    #file_str.write("v {} {} {} {}\n".format(p[0], p[1], p[2], p[3]))

                    uv = np.array((u_base + j / degu * (u_end - u_base), v_base + i / degv * (v_end - v_base)))
                    uv -= np.array(uvbbox[0])
                    uv[0] /= (uvbbox[1][0] - uvbbox[0][0])
                    uv[1] /= (uvbbox[1][1] - uvbbox[0][1])
                    
                    #file_str.write("vt {} {}\n".format(uv[0], uv[1]))

                    bcpts[i][j] = np.array([p[0], p[1], p[2], p[3], uv[0], uv[1]])

                    # print(xyz,w,uv)

            # convert bezier patch to bicubic
            vpiece = 1 if degv <= 3 else 3
            upiece = 1 if degu <= 3 else 3

            _bcpts = convert2Cubic([bcpts[:,i] for i in range(degu+1)])
            bcpts = np.array(_bcpts)
            _bcpts = convert2Cubic([bcpts[:,i] for i in range(degv+1)])
            bcpts = np.array(_bcpts)

            for i in range(vpiece):
                for j in range(upiece):
                    total_bicubic_piece += 1
                    for k in range(4):
                        for l in range(4):
                            p = bcpts[i*3 + k][j*3 + l]
                            
                            file_str.write("v {} {} {} {}\n".format(p[0], p[1], p[2], p[3]))
                            file_str.write("vt {} {}\n".format(p[4], p[5]))
                    
            
                    
    idx = 1
    # for piece_v in range((cpts.shape[0] - 1)//degv):
    #     for piece_u in range((cpts.shape[1] - 1)//degu):
    #         file_str.write("p {} {}\n".format(degu, degv))
    #         for i in range((degu + 1)*(degv + 1)):
    #             file_str.write("{}/{}\n".format(idx, idx))
    for i in  range(total_bicubic_piece):
        file_str.write("p {} {}\n".format(3, 3))
        for i in range(16):
            file_str.write("{}/{}\n".format(idx, idx))
            idx+=1

    return file_str.getvalue()

def convert2Cubic(cpts):
    if len(cpts) == 2:
        return degRaiseFrom1(cpts)
    elif len(cpts) == 3:
        return degRaiseFrom2(cpts)
    elif len(cpts) > 4:
        return degLower2Cubic(cpts)
    else:
        return cpts

def degRaiseFrom1(cpts):
    new_cpts = [None] * 4
    new_cpts[0] = cpts[0]
    new_cpts[1] = (cpts[0] * 2.0 + cpts[1])/ 3.0
    new_cpts[2] = (cpts[0] + 2.0 * cpts[1])/ 3.0
    new_cpts[3] = cpts[1]
    return new_cpts
    
def degRaiseFrom2(cpts):
    new_cpts = [None] * 4
    new_cpts[0] = cpts[0]
    new_cpts[1] = (cpts[0] + 2.0 * cpts[1])/ 3.0
    new_cpts[2] = (cpts[2] + 2.0 * cpts[1])/ 3.0
    new_cpts[3] = cpts[2]
    return new_cpts

def degLower2Cubic(cpts):
    """ degree lower a high degree bezier funtion to 3 cubic function
    """

    deg = len(cpts) - 1
    assert(deg >= 1)

    start_points = [None]*10

    # get the first 3 control points by 3 times de Castljue
    start_points[0] = cpts[0]
    start_points[1] = cpts[0] * 2/3 + cpts[1] * 1/3
    start_points[2] = cpts[0] * 4/9 + cpts[1] * 4/9 + cpts[2] * 1/9
    
    # get the last 3 control points
    start_points[9] = cpts[-1]
    start_points[8] = cpts[-1] * 2/3 + cpts[-2] * 1/3
    start_points[7] = cpts[-1] * 4/9 + cpts[-2] * 4/9 + cpts[-3] * 1/9

    # calculate the middle 4 points
    start_points[3] = 1/3 * start_points[7] + 7/6 * start_points[2] - 1/6 * start_points[8] - 1/3 * start_points[1]
    start_points[4] = 2/3 * start_points[7] + 4/3 * start_points[2] - 1/3 * start_points[8] - 2/3 * start_points[1]
    start_points[5] = 2/3 * start_points[2] - 1/3 * start_points[1] + 4/3 * start_points[7] - 2/3 * start_points[8]
    start_points[6] = 1/3 * start_points[2] - 1/6 * start_points[1] + 7/6 * start_points[7] - 1/3 * start_points[8]

    return start_points
    

def exportNURBSCurve(curve, bbox=[(0,0), (1,1)]):
    file_str = io.StringIO()
    
    degree = curve.M
    assert(degree <= 3)
    assert(curve.prop3 == 1)

    # print(curve)

    # first convert nurbs to bezier
    # TODO: assume the knot sequence has knot multiplicty degree + 1 at both end
    b_knot, b_cpts = nurbs2bezier(curve.T[1:-1], [np.array(p) for p in curve.control_points], degree)

    for p in b_cpts:
            p[0:2] -= np.array(bbox[0])
            p[0] /= bbox[1][0] - bbox[0][0]
            p[1] /= bbox[1][1] - bbox[0][1]
        

    file_str.write("M {},{}".format(b_cpts[0][0], b_cpts[0][1]))

    first = True
    
    for i in range((len(b_cpts) - 1) // degree):
        print("-----")

        cpts = []
        for j in range(degree + 1):
            newp = b_cpts[i*degree + j]
            cpts.append(newp)

        if degree == 1:
            file_str.write(" L {},{}".format(cpts[1][0], cpts[1][1]))
            continue
        
        if degree == 2:
            cpts = degRaiseFrom2(cpts)

        for i,p in enumerate(cpts):
            if first:
                first = False
                file_str.write(" C ")
            elif i != 0:
                file_str.write(" {},{}".format(p[0], p[1]))

    return  file_str.getvalue()
        

# # test degLoer2Cubic
# cpts = [0,1,2,3,4,5]

# print(degLower2Cubic(cpts))

# sys.exit(0)
    
idx = 0
for entity in entity_list:
    # if idx >= 1:
    #      break;
    if entity.d['entity_type_number'] == 143:
        print(entity)
        # print("sptr: ", entity.SPTR)
        surface = entity_list[pointer_dict[entity.SPTR]]

        # only support nurbs surface for now
        if surface.d['entity_type_number'] != 128:
            continue
            
        if surface.M1 == 3 and surface.M2 == 3:
            print("export!")
            print(surface)
            idx += 1
            # if idx != 61:
            #      continue
            with open("trim_surface_{}.txt".format(idx), "w") as f:
                # print(surface)
                f.write(exportNURBSSurface(surface))
            bbox = [(surface.U0, surface.V0), (surface.U1, surface.V1)]
            
            with open("trim_curve_{}.txt".format(idx), "w") as f:
                for bptr in entity.BDPT:
                    boundary = entity_list[pointer_dict[bptr]]
                    print(boundary)
                    for crvpt, sense, pscpts in boundary.PSCPT:
                        for pscpt in pscpts:
                            pcurve = entity_list[pointer_dict[pscpt]]
                            # print(exportNURBSCurve(pcurve))
                            f.write(exportNURBSCurve(pcurve, bbox))
                            if pcurve.prop2 == 1:
                                f.write(" z")
                            f.write(" ")
        
    elif entity.d['entity_type_number'] == 144:
        # print(entity)
        surface = entity_list[pointer_dict[entity.PTS]]
        if surface.d['entity_type_number'] != 128:
            continue
            
        
        if True:
            print("export!")
            #print(surface)
            
            idx += 1
            # if idx != 173:
            #      continue
            with open("trim_surface_{}.txt".format(idx), "w") as f:
                print(surface)
                f.write(exportNURBSSurface(surface))
            bbox = [(surface.U0, surface.V0), (surface.U1, surface.V1)]
            
            with open("trim_curve_{}.txt".format(idx), "w") as f:
                CPTS = []

                flipped = False

                if entity.N1 == 0:
                    flipped = True
                else:
                    CPTS.append(entity.PTO)

                assert(flipped == False)

                CPTS.extend(entity.PTI)
                
                for bptr in CPTS:
                    pcurve = entity_list[pointer_dict[bptr]]
                    curve = entity_list[pointer_dict[pcurve.BPTR]]

                    assert(curve.d['entity_type_number'] == 102)
                    print(curve)
                    for de in curve.DE:
                        e = entity_list[pointer_dict[de]]
                        assert(e.d['entity_type_number'] == 126)
                        print(e)
                    
                        f.write(exportNURBSCurve(e, bbox))
                        if e.prop2 == 1:
                             f.write(" z")
                        f.write(" ")

        else:
            pass
            #print("not support!")
            #print(surface.M1, surface.M2)
        
                        
    #elif entity.d['entity_type_number'] == 128:
    #     print(entity)
