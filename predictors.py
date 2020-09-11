import numpy as np


def rnn_predict(x0, v0, predictor, steps0, steps, mr, dt, t0=0, option=0):
    
    ## option = 0 if model predicts v - vm; option = 1 if model predicts v - u
    assert option == 0 or option == 1, "option not recognized!"

    ## run through spinup stage
    xs, _, us, Dudts, _ = fullmr_predict(x0, v0, steps0, mr, dt, t0)
    feats = np.concatenate((us, Dudts), axis=-1)
    feats = (feats - predictor.input_mean)/predictor.input_std
    
    if feats.ndim < 3:
        feats = np.expand_dims(feats, 0)
    else:
        feats = np.swapaxes(feats, 0, -2)

    x0 = xs[-1,...,:]
    u0 = us[-1,...,:]
    t = t0 + dt*steps0    #### t0 should be set before the spinup stage
    outputs = predictor.spinup(feats)
    vd, lstm_states = np.squeeze(outputs[0]), outputs[1:]
    if option == 0:
        v0 = mr.order1_manifold_vel(x0, t) + vd
    elif option == 1:
        v0 = u0 + vd
    
    # quantities to record
    x, v, u, DuDt, vm = [], [], [], [], []
    vt = []
    
    ## main time-march loop
    for _ in np.arange(steps):
        
        x1 = x0 + dt*v0
        u1 = mr.flow.vel(x1, t+dt)
        DuDt1 = mr.flow.DuDt(x1, t+dt)
        
        feats1 = (np.concatenate((u1, DuDt1), axis=-1) - predictor.input_mean)/predictor.input_std
        feats1 = np.expand_dims(feats1, axis=-2)    # expand once for time
        if feats1.ndim < 3:
            feats1 = np.expand_dims(feats1, 0)
        
        outputs = predictor.step([feats1]+lstm_states)
        vd, lstm_states = np.squeeze(outputs[0]), outputs[1:]
        vm1 = mr.order1_manifold_vel(x1, t+dt)
        
        if option == 0:
            v1 = vm1 + vd
        elif option == 1:
            v1 = u1 + vd

        # calculate true particle velocity following the MODEL trajectory
        _, vt1 = mr.RK4step(x0, v0, dt, t)

        x.append(x1)
        u.append(u1)
        DuDt.append(DuDt1)
        v.append(v1)
        vm.append(vm1)
        vt.append(vt1)
        
        x0, v0 = x1, v1
        t += dt
        
    return np.array(x), np.array(v), np.array(u), np.array(DuDt), np.array(vm), np.array(vt)
    

def fullmr_predict(x0, v0, steps, mr, dt, t0=0):
    
    t = t0
    x, v, u, DuDt, vm = [], [], [], [], []    # quantities to record
    
    for _ in np.arange(steps):
    
        x1, v1 = mr.RK4step(x0, v0, dt, t)
        u1 = mr.flow.vel(x1, t+dt)
        DuDt1 = mr.flow.DuDt(x1, t+dt)
        vm1 = u1 + mr.ep*(3/2*mr.R - 1)*DuDt1
        
        # record
        x.append(x1)
        v.append(v1)
        u.append(u1)
        DuDt.append(DuDt1)
        vm.append(vm1)
        
        # update
        x0, v0 = x1, v1
        t += dt
        
    return np.array(x), np.array(v), np.array(u), np.array(DuDt), np.array(vm)    # 2D arrays with dim [steps, 2]


def slowmanifold_predict(x0, steps, mr, dt, t0=0):
    # t0 should be set after the spinup stage
    t = t0
    vm0 = mr.order1_manifold_vel(x0, t)
    x, u, DuDt, vm = [], [], [], []    # quantities to record
    
    for _ in np.arange(steps):
    
        x1 = x0 + dt*vm0
        u1 = mr.flow.vel(x1, t+dt)
        DuDt1 = mr.flow.DuDt(x1, t+dt)
        vm1 = u1 + mr.ep*(3/2*mr.R - 1)*DuDt1
        
        # record
        x.append(x1)
        u.append(u1)
        DuDt.append(DuDt1)
        vm.append(vm1)
        
        # update
        x0, vm0 = x1, vm1
        t += dt
        
    return np.array(x), np.array(u), np.array(DuDt), np.array(vm)    # 2D arrays with dim [steps, 2]
