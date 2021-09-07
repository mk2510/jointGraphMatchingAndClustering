import numpy as np 


def groupLs(ls):
    return_ls = []
    index = 0
    size = 1
    while index < len(ls):
        groupFinal = False
        size = 1

        while not groupFinal:
            ref_group = np.array([ls[index]] * size)
            comp_group = ls[index:index + size]

            result = np.allclose(ref_group, comp_group, atol=0.1)

            if result:
                size += 1
                if index + size > len(ls):
                    size -= 1
                    groupFinal = True
            else:
                groupFinal = True
                size -= 1
        
        return_ls.append(size)
        index += size
    return np.array(return_ls)
