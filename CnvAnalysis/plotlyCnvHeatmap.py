import chart_studio.plotly as py
import plotly
import plotly.graph_objs as go
from sys import stdin,argv
binsize = int(argv[1])

chromosomes=['1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20','21','22','X']
#use hg19
chromsize=[249250621,243199373,198022430,191154276,180915260,
171115067,159138663,146364022,141213431,135534747,
135006516,133851895,115169878,107349540,102531392,
90354753,81195210,78077248,59128983,63025520,
48129895,51304566,155270560]

x_axis = []
x_ticks=[]
x_ticktext=[]
for i in range(len(chromosomes)):
    x_ticks=x_ticks+['chr'+chromosomes[i]+'-1']
    x_ticktext=x_ticktext+['chr'+chromosomes[i]]
    #x_ticktext = x_ticktext + ['']
    x_axis = x_axis+['chr'+chromosomes[i]+'-'+str(j+1) for j in range(int(chromsize[i]/binsize)+1)]
z =  []
samples=[]
for line in stdin:
    line=line.strip()
    a = line.split('\t')
    # if sample name is at a[-1] run the following codes
    #b = [float(i) for i in a[0:len(a)-1]]
    #z.append(b)
    #samples.append(a[-1])

    # if sample names defined separately run the following codes
    b = [float(i) for i in a]
    z.append(b)


samples = ['Ln7','Ln9','Ln1','BrM','BrP',
             'Ln11','Ly2','Ln3',
             'Bo3','Ln10','Bo1','Ln8','Lv3','Ln5','Bo2','Bn2','Bn1','Bn3','Bn4','Ln2',
             'Ly1','Ln6',
             'Kd1','Ln4','Lv4','Lv2','Lv1','Pa1']

#data = [go.Heatmap(z=z, x=x_axis,y=samples, colorscale=[[1,'rgb(255,0,0)'],[0,'rgb(0,0,255)'],[0.5,'rgb(255,255,255)']])]
data = [go.Heatmap(z=z, x=x_axis,y=samples, colorscale='Picnic',zauto=False, zmin=-2,zmax=2)]
layout = go.Layout(
    title='CNV',
    xaxis = dict(ticks='inside',tickmode='array',tickvals=x_ticks,ticktext=x_ticktext,ticklen=1037),
    yaxis = dict(ticks='')
)


#print(len(x_axis))
#print(len(z[0]))
#print('x_axis', x_axis)
#print(z)


fig = go.Figure(data=data, layout=layout)
plotly.offline.plot(fig, filename='cnv-heatmap.html')
