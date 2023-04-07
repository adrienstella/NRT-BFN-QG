import numpy as np

def h2uv(ssh,dy,dx,g,f, ubc=None,vbc=None):
    """ SSH to U,V

    Args:
        ssh (2D array): SSH field.
        dy (2D array): Height of grid points
        dx (2D array): Width of grid points
        f (2D array): Coriolis parameter
        g (scalar): Acceleration of gravity (m.s-2)

    Returns:
        u (2D array): Zonal velocity  
        v (2D array): Meridional velocity
    """

    h = ssh.values

    ny = len(ssh.lat)
    nx = len(ssh.lon)

    u = np.zeros((ny,nx))
    v = np.zeros((ny,nx))

    u[1:-1,1:] = - g/f[1:-1,1:]*(h[2:,:-1]+h[2:,1:]-h[:-2,1:]-h[:-2,:-1])/(4*dy[1:-1,1:])
    v[1:,1:-1] = + g/f[1:,1:-1]*(h[1:,2:]+h[:-1,2:]-h[:-1,:-2]-h[1:,:-2])/(4*dx[1:,1:-1])

    # Condition de Neumann ( du/dn = 0 )
    u[:,0] = u[:,1]
    u[0,:] = u[1,:]
    u[-1,:] = u[-2,:]
    
    v[0,:] = v[1,:]
    v[:,0] = v[:,1]
    v[:,-1] = v[:,-2]
    
    return u,v


def h2rv(ssh,dy,dx,g,f):
    """ SSH to Q

    Args:
        h (2D array): SSH field.
        dy (2D array): Height of grid points
        dx (2D array): Width of grid points
        f (2D array): Coriolis parameter
        g (scalar): Acceleration of gravity (m.s-2)

    Returns:
        xi_norm: Normalized relative vorticity field  
    """
        
    h = ssh.values

    ny = len(ssh.lat)
    nx = len(ssh.lon)

    xi_norm = np.zeros((ny,nx))

    # Normalized relative vorticity
    xi_norm[1:-1,1:-1] = (1/f[1:-1,1:-1]) * g/f[1:-1,1:-1]*\
        ((h[2:,1:-1]+h[:-2,1:-1]-2*h[1:-1,1:-1])/dy[1:-1,1:-1]**2 + (h[1:-1,2:]+h[1:-1,:-2]-2*h[1:-1,1:-1])/dx[1:-1,1:-1]**2)

    # Condition de Neumann ( dq/dn = 0 )
    xi_norm[:,0] = xi_norm[:,1]
    xi_norm[:,-1] = xi_norm[:,-2]
    xi_norm[0,:] = xi_norm[1,:]
    xi_norm[-1,:] = xi_norm[-2,:]

    return xi_norm


def h2pv(ssh,dy,dx,g,f, c=None):
    """ SSH to Q

    Args:
        h (2D array): SSH field.
        c (2D array): Phase speed of first baroclinic radius 
        dy (2D array): Height of grid points
        dx (2D array): Width of grid points
        f (2D array): Coriolis parameter
        g (scalar): Acceleration of gravity (m.s-2)

    Returns:
        q: Potential Vorticity field  
    """
        
    h = ssh.values

    ny = len(ssh.lat)
    nx = len(ssh.lon)

    q = np.zeros((ny,nx))

    # Potential vorticity
    q[1:-1,1:-1] = g/f[1:-1,1:-1]*\
        ((h[2:,1:-1]+h[:-2,1:-1]-2*h[1:-1,1:-1])/dy[1:-1,1:-1]**2 +\
            (h[1:-1,2:]+h[1:-1,:-2]-2*h[1:-1,1:-1])/dx[1:-1,1:-1]**2) -\
            g*f[1:-1,1:-1]/(c[1:-1,1:-1]**2) *h[1:-1,1:-1]

    # Condition de Neumann ( dq/dn = 0 )
    q[:,0] = q[:,1]
    q[:,-1] = q[:,-2]
    q[0,:] = q[1,:]
    q[-1,:] = q[-2,:]

    return q