import numpy as np
import h5py



class GustoDataset(object):

    def __init__(self, filename):
        self.h5f = h5py.File(filename, 'r')

    def close(self):
        self.h5f.close()

    def get_vertex_positions(self):
        """ Return the array: verts[row][col][xyz] """
        verts = [ ]
        for row in self.h5f['rows'].itervalues():
            x1 = row['verts']['x1'][:]
            x2 = row['verts']['x2'][:]
            x3 = row['verts']['x3'][:]
            coords = np.vstack([x1, x2, x3]).T
            verts.append(coords)
        return verts

    def get_cell_variable(self, key):
        cells = [ ]
        for row in self.h5f['rows'].itervalues():
            var = row['cells'][key][:]
            cells += var.flat
        return np.array(cells)

    def get_vert_variable(self, key, flat=True):
        verts = [ ]
        for row in self.h5f['rows'].itervalues():
            var = row['verts'][key][:]
            verts.append(var)
        if flat:
            return np.array(verts).flatten()
        else:
            return verts

    def get_face_segments(self):
        verts = self.get_vertex_positions()
        faces = self.h5f['faces'][:]
        segments = [ ]
        for face in faces:
            v0 = face[0]
            v1 = face[1]
            row_index0, col_index0 = v0
            row_index1, col_index1 = v1
            seg = [verts[row_index0][col_index0],
                   verts[row_index1][col_index1]]
            segments.append(seg)
        return np.array(segments)

    def get_cell_polygons(self):
        polygons = [ ]
        for row in self.h5f['rows'].itervalues():
            cells = row['cells']
            v0_x1 = cells['v0.x1'][:]
            v0_x2 = cells['v0.x2'][:]
            v0_x3 = cells['v0.x3'][:]
            v1_x1 = cells['v1.x1'][:]
            v1_x2 = cells['v1.x2'][:]
            v1_x3 = cells['v1.x3'][:]
            v2_x1 = cells['v2.x1'][:]
            v2_x2 = cells['v2.x2'][:]
            v2_x3 = cells['v2.x3'][:]
            v3_x1 = cells['v3.x1'][:]
            v3_x2 = cells['v3.x2'][:]
            v3_x3 = cells['v3.x3'][:]
            v0 = zip(v0_x1, v0_x3)
            v1 = zip(v1_x1, v1_x3)
            v2 = zip(v2_x1, v2_x3)
            v3 = zip(v3_x1, v3_x3)
            polygons += list(zip(v0, v1, v3, v2))
        return np.array(polygons)
