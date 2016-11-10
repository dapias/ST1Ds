"""
Given a 3x3 matrix  orthogonalize its columns with the Gram-Schmidt procedure
"""
function gramschmidt(u::Matrix{Float64})
     w = eye(3)
     w[:,1] = u[:,1];
     v1 = w[:,1]/norm(w[:,1])
     w[:,2] = u[:,2] - dot(u[:,2],v1)*v1;
     v2 = w[:,2]/norm(w[:,2]);
     w[:,3] = (u[:,3] - dot(u[:,3],v2)*v2 - dot(u[:,3],v1)*v1)
    
    return w
end
