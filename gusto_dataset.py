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
        for row in self.h5f['rows']:
            x1 = self.h5f['rows'][row]['verts']['x1'][:]
            x2 = self.h5f['rows'][row]['verts']['x2'][:]
            x3 = self.h5f['rows'][row]['verts']['x3'][:]
            coords = np.vstack([x1, x2, x3]).T
            verts.append(coords)
        return verts

    def get_cell_variable(self, key):
        cells = [ ]
        for row in self.h5f['rows']:
            var = self.h5f['rows'][row]['cells'][key][:]
            cells.append(var)
        return cells

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
