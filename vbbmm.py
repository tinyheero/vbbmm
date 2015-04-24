#!/usr/bin/env python

from __future__ import division

import csv
import numpy as np
import sys

from pyvmp.inference_engine import InferenceEngine
from pyvmp.nodes.vector import Dirichlet, Discrete, Multinomial
from pyvmp.nodes.mixture import MixtureDistribution
from pyvmp.graph import ModelGraph

def main(args):
    data, positions = load_data(args.in_file, args.num_samples, min_depth=args.min_depth)
    
    print 'Clustering {0} data points'.format(len(data))
    
    mix_weight_priors = [args.mix_weight_priors] * args.num_components
    
    component_priors = [[1, 1], ] * args.num_components

    results = binomial_mixture_model(data, component_priors, mix_weight_priors, args.num_samples, max_iters=args.max_iters)
    
    labels = []
    
    for x in results:
        labels.append(np.argmax(x[0].mix_weights))
    
    labels = np.array(labels)

    with open(args.out_file, 'w') as out_fh:
        writer = csv.DictWriter(out_fh, ['chr', 'pos', 'cluster'], delimiter='\t')
        
        writer.writeheader()
    
        for (chrom, coord), l in zip(positions, labels):
            writer.writerow({'chr' : chrom, 'pos' : coord, 'cluster' : l})

def load_data(file_name, num_samples, min_depth=0):
    data = []
    
    positions = []
    
    reader = csv.DictReader(open(file_name), delimiter='\t')
    
    cols = []
    
    for i in range(1, num_samples + 1):
        cols.append('T{0}_TR'.format(i))
        
        cols.append('T{0}_TA'.format(i))
    
    print 'Analysing data from columns {0}'.format(','.join(cols))
    
    for row in reader:
        data_row = []
 
        for col in cols:
            data_row.append(int(row[col]))
        
        passed_depth_threshold = True
        
        for i in range(num_samples):
            sample_data = data_row[2 * i : 2 * (i + 1)] 

            if sum(sample_data) < min_depth:
                passed_depth_threshold = False
        
        if passed_depth_threshold:
            data.append(data_row)
            
            positions.append((row['chr'], row['pos']))

    return np.array(data), positions

def binomial_mixture_model(data, density_priors, mix_weight_priors, num_components, max_iters=100):
        
    print 'Initialising model'
    
    # A container for the nodes in the model.
    model = ModelGraph()
    
    mix_weight_priors = load_mix_weight_priors(mix_weight_priors, model)
    
    component_priors = []
    
    for _ in range(num_components):
        component_priors.append(load_component_priors(density_priors, model))
    
    x = []
    
    # Add data points
    for row in data:
        index = Discrete({'pi' : mix_weight_priors})
        
        model.add_node(index)
        
        x_i_nodes = []
        
        for i, comp_prior in enumerate(component_priors):
            parents = {}
            
            parents['index'] = index
            
            parents['components'] = comp_prior
        
            node = MixtureDistribution(Multinomial, parents, value=row[2 * i:  2 * (i + 1)])    
        
            model.add_node(node)
            
            x_i_nodes.append(node)
        
        x.append(x_i_nodes)

    # Helps convergence if we update the Dirichlet and indicator node first.
    mix_weight_priors.update()
    
    for node in model.nodes:
        if isinstance(node, Discrete):
            node.update()            
    
        
    print 'Fitting model'
    
    # Create the inference engine and run.
    ie = InferenceEngine(model, verbose=True)
    
    ie.run(max_iters=max_iters)

    print 'Fitting finished'    

    print "#" * 100
    print mix_weight_priors.e_x

    print "#" * 100    
    for comp_prior in component_priors:
        for comp_parents in comp_prior:
            print comp_parents['pi'].e_x
    
    return x
        
def load_component_priors(priors, model):
    component_priors = []

    for _, alpha in enumerate(priors):
        prior = Dirichlet({'alpha' : alpha})
        
        model.add_node(prior)
        
        component_priors.append({'pi' : prior})
        
    return component_priors

def load_mix_weight_priors(mix_weight_priors, model):
    # Mix-weight prior
    mix_weight_priors = Dirichlet({'alpha' : mix_weight_priors})    
    
    model.add_node(mix_weight_priors)
    
    return mix_weight_priors

if __name__ == '__main__':
    import argparse
    
    parser = argparse.ArgumentParser()
    
    parser.add_argument('--in_file', required=True)
    
    parser.add_argument('--out_file', required=True)
    
    parser.add_argument('--max_iters', default=100, type=int,
                        help='''Maximum number of VB iterations to do. Default is 100.''')    
    
    parser.add_argument('--min_depth', default=0, type=int,
                        help='''Minimum depth in any sample required to include position for analysis.''')
    
    parser.add_argument('--mix_weight_priors', default=1e-6, type=float,
                        help='''Dirichlet prior on component mix-weights. Default is 1e-6.''')
    
    parser.add_argument('--num_components', default=10, type=int,
                        help='''Number of mixture components to use. Default is 10.''')
    
    parser.add_argument('--num_samples', default=2, type=int,
                        help='''Number of samples to analyses. Default is 2.''')
    
    args = parser.parse_args()
    
    main(args)
