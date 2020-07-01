import numpy as np


def cwbrain(plotly=False):
    """Rain color bar for Central Weather Bureau in Taiwan"""
    colors = np.array([
     [1.000,1.000,1.000], [0.000,0.000,0.000], \
     [1.000,1.000,1.000], [0.608,1.000,1.000], \
     [0.000,0.812,1.000], [0.039,0.596,1.000], [0.039,0.396,1.000], \
     [0.188,0.600,0.039], [0.196,1.000,0.000], \
     [0.973,1.000,0.000], [1.000,0.796,0.000], [1.000,0.603,0.000], \
     [0.980,0.012,0.000], [0.800,0.000,0.012], [0.627,0.000,0.000], \
     [0.596,0.000,0.604], [0.765,0.016,0.800], \
     [0.973,0.020,0.953], [0.996,0.796,1.000], \
    ], 'f')
    if plotly:
       colors = np.round(colors*255).astype(int)
       plotly_type = []
       i = 0.0
       colorbin = 1.0/float(len(colors)-2-1)
       for j,color in enumerate(colors):
           if j > 1:
              #print("--{}".format(str(color)[1:12]))
              #print(i) # 0 ~ 1
              plotly_type.append([i,'rgb({}, {}, {})'\
                     .format(color[0],color[1],color[2])])
              i += colorbin
       #print("--{}".format(plotly_type))
       colors = plotly_type
    return colors


if __name__ == '__main__':
   cwbrain()
