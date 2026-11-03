from pathlib import Path
import sys
import numpy as np
import pytest
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'code'))
from registration_geometry import jacobian_from_displacement,symmetric_surface_distances

def test_oblique_physical_affine_jacobian():
    theta=.4;rot=np.array([[np.cos(theta),-np.sin(theta),0],[np.sin(theta),np.cos(theta),0],[0,0,1]])
    grid=rot@np.diag([2,3,4]);mapping=np.array([[1.1,.2,.1],[0,.9,.05],[0,0,1.2]])
    z,y,x=np.indices((6,7,8));p=np.stack([x,y,z],axis=-1)@grid.T
    u=p@(mapping-np.eye(3)).T
    actual=jacobian_from_displacement(u,grid)
    np.testing.assert_allclose(actual,np.linalg.det(mapping),atol=2e-7,rtol=2e-7)

def test_detects_reflection_and_identity():
    z,y,x=np.indices((5,5,5));u=np.zeros((5,5,5,3));assert np.all(jacobian_from_displacement(u,np.eye(3))==1)
    u[...,0]=-2*x;np.testing.assert_allclose(jacobian_from_displacement(u,np.eye(3)),-1)

def test_surface_distance_is_in_mm_and_identical_masks_zero():
    a=np.zeros((9,9,9),bool);a[4,4,3]=1;b=np.zeros_like(a);b[4,4,4]=1;domain=np.ones_like(a)
    r=symmetric_surface_distances(a,b,[3,2,2],domain);assert r['mean_symmetric_mm']==pytest.approx(2)
    assert symmetric_surface_distances(a,a,[3,2,2],domain)['HD95_mm']==0
