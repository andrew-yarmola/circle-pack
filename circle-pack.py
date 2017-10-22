#!/usr/bin/env python3

from PyQt5.QtCore import (Qt, QPointF, QLineF, QSize)
from PyQt5.QtGui import (QPainter, QPainterPath, QValidator, QColor)
from PyQt5.QtWidgets import (QApplication, QMainWindow, QWidget,
                             QPushButton, QCheckBox, QLabel,
                             QLineEdit, QGridLayout)

import numpy as np
from numpy import exp, pi, sqrt, sin, cos, log, arccosh
from collections import defaultdict

def circle(T) :
    z1,z2,z3 = T
    ax, ay = z1.real, z1.imag
    bx, by = z2.real, z2.imag
    cx, cy = z3.real, z3.imag
    a = bx - ax
    b = by - ay
    c = cx - ax
    d = cy - ay
    e = a * (ax + bx) + b*(ay + by)
    f = c * (ax + cx) + d*(ay + cy)
    g = 2 * (a * (cy - by) - b * (cx - bx))
    if g == 0. :
        return (z1, z2)
    else :
        center = (d*e - b*f)/g + 1j *(a*f - c*e)/g
        radius = sqrt((ax - center.real)**2 + (ay - center.imag)**2)
        return (center, radius)

def shear_graft(s,x,y,z) :
    # Multiplicative shear s and graft pi/2 accross (x,y) edge.
    # Returns 4th point, opposite z.
    return ((1j + s)*x*y - 1j*x*z - s*y*z)/(s*x + 1j*y - (1j + s)*z)

def shear(s,x,y,z) :
    # Multiplicative shear s accross (x,y) edge.
    # Returns 4th point, opposite z.
    return ((1 + s)*x*y - (x + s*y)*z)/(s*x + y - (1 + s)*z)

def L(T,s) :
    return (shear(s,T[0],T[1],T[2]),T[1],T[0])

def R(T,s) :
    return (shear(s,T[2],T[0],T[1]),T[0],T[2])

def LGr(T,s) :
    return (shear_graft(s,T[0],T[1],T[2]),T[1],T[0])

def RGr(T,s) :
    return (shear_graft(s,T[2],T[0],T[1]),T[0],T[2])

def rot_L(T) :
    return (T[2],T[0],T[1])

def rot_R(T) :
    return (T[1], T[2], T[0])

def shear_up(s,t) :
    return [t, (-1+2*(s-t)*t)/(s-2*t+2*s*t**2), (s-2*t+2*s*t**2)/(-1+2*(s-t)*t), t, 1]

def shear_right(s,t) :
    return [1/t, s, 1/s, 1/t, 1]

def trace_one(s,t) :
    return (-1j - s + 2*t - 2*(s - 1j)*t**2)/(s*t)

def trace_two(s,t) :
    return (-1j - s + 2*1j*s*t - 4*1j*s*t**3 - 4*(s - 1j)*t**4)/(t*(s - 2*t + 2*s*t**2))

def conj_param_squared(s,t) :
    return (s - 2*1j*s*t**2 + 1j*(1 + (1 + 1j)*t)**2)/(1 - 1j*s - (2 - 2*1j)*t + 2*(s - 1j)*t**2)

# We can assume from the packing combinatorics that |Im(c)| <= Pi and
# |Im(c tau) | <= Pi because the big circle is self tangent so its translates
# can't go too far away. The equations below are adapted to that the imaginary
# satify the restrictions.

def c_param(s,t) :
    return 2*arccosh(trace_one(s,t)/2)

# unused
def c_tau_param(s,t) :
    return 2*arccosh(trace_two(s,t)/2)

def tau_param(s,t) :
    # slightly different from Mathematica. Possibly an difference in arccosh implementation
    if t > 1/sqrt(2) :
        return -arccosh(trace_two(s,t)/2)/arccosh(trace_one(s,t)/2)
    else :
        return arccosh(trace_two(s,t)/2)/arccosh(trace_one(s,t)/2)

def fern(W,H, sR, sU) :
    # 2W and 2H will be the total grid of the 4 pt circle
    # s_ are shear coords seqeunces for right and up generators.
    # we start with the first key triangles, the first one will always be the top one.
    # half coordinates mean intermediate triangles, integer are pt circles
    # always "top, bottom" triangle
    all_T = defaultdict(list)
    all_T[(0,0)] = [(exp(3*pi*1j/4), exp(-3*pi*1j/4), exp(pi*1j/4)),
                    (exp(-pi*1j/4), exp(pi*1j/4),exp(-3*pi*1j/4))]
    paths = {('vertical',1) : [RGr,RGr,LGr,RGr,L], ('vertical',-1) : [RGr,LGr,RGr,LGr,L],
            ('horizontal',1) : [LGr,LGr,RGr,LGr,R], ('horizontal',-1) : [LGr,RGr,LGr,RGr,R]}
    s_vertical = {1 : sU, -1 : sU[-2::-1] + [sU[-1]], 1/2 : sU[0], -1/2 : sU[-2]}
    s_horizontal = {1 : sR, -1 : sR[-2::-1] + [sU[-1]]}
    
    for sgn in [-1,1] :
        path = paths[('vertical',sgn)]
        shears = s_vertical[sgn]
        for m in range(H) :
            T = all_T[(0, m*sgn)][(1-sgn)//2]
            new = []
            for k in range(len(path)) :
                T = path[k](T, shears[k])
                if k != 1 :
                    new.append(T)
            all_T[(0,(m + 1/2)*sgn)] = new[:3]
            all_T[(0,(m + 1)*sgn)] = [rot_L(new[-2]), new[-1]][::-sgn]

    for m in range(-H,H+1) :
        for sgn in [-1,1] :
            path = paths[('horizontal',sgn)]
            shears = s_horizontal[sgn]
            for n in range(W) :
                T = all_T[(n*sgn,m)][(1+sgn)//2]
                new = []
                for f,s in zip(path,shears) :
                    T = f(T,s)
                    new.append(T)
                all_T[((n + 1/2)*sgn, m)] = new[:3]
                all_T[((n + 1)*sgn, m)] = [rot_R(new[-2]), new[-1]][::sgn]
                for ud in [-1,1] :
                    T = all_T[((n + 1)*sgn, m)][(1-ud)//2]
                    all_T[((n + 1)*sgn, m + ud/2)].append(RGr(T, s_vertical[ud/2]))

    return all_T

def all_circles(all_T, include_dual, x = None) :
    all_circles = {'small' : [], 'large' : [], 'dual' : []}
    for p, v in all_T.items() :
        x_whole = (p[0] - int(p[0]) == 0)
        y_whole = (p[1] - int(p[1]) == 0)
        if x is not None :
            v = tuple(map(lambda z : (x+z)/(x-z), v))
        if x_whole and y_whole :
            all_circles['small'].append(circle(v[0]))
        else :
            if y_whole and not x_whole :
                all_circles['large'].append(circle(v[1]))
                if include_dual :
                    all_circles['dual'].extend((circle(v[0]),circle(v[2])))
            elif include_dual :
                all_circles['dual'].extend(map(circle,v))

    return all_circles

from scipy.linalg import solve_banded

def bezier_control_points(data) :
    # Assume data is a numpy arry of shape (k,2)
    count = len(data) - 1
    if count < 1 :
        raise Exception("Requires at least two points")
    
    # Our equations are  :
    # 2 P1_0 + P1_1 = D_0 + 2 D_1
    # P1_{i-1} + 4 P1_i + 1 P1_{i+1} = 4 D_i + 2 D_{i+1}
    # 2 P1_{n-2} + 7 P1_{n-1} = 8 D_{n-1} + D_n
    
    diag = np.tile([1,4,1],(count,1))
    diag[0] = [0,2,1]
    diag[-1] = [1,7,0]
    diag[-2,2] = 2
    
    rhs = 4 * data[:count] + 2 * data[1:]
    rhs[0] -= 3 * data[0]
    rhs[-1] += 4 * data[count-1] - data[count]
    
    P1 = solve_banded((1,1), diag.transpose(), rhs,
                      overwrite_ab = True, overwrite_b = True,
                      check_finite = False)
    # Equations for P2
    # P2_i = 2 D_{i+1} - P1_{i+1}
    # P_{n-1} = (D_n + P1_{n-1})/2
    P2 = np.zeros((count,2))
    P2[:-1] = 2 * data[1:-1] - P1[1:]
    P2[-1] = (data[-1] + P1[-1])/2
                      
    return (P1, P2)

def _get_point(evt):
    return (evt.pos().x(),evt.pos().y())

from itertools import starmap

class ParamSpace(QWidget) :

    drag_dist = 4
    max_t = 2.5
    scale = (400.,195.)
    shift = (20.,490.)
    default_param = (sqrt(2.) - 1., 1./sqrt(2))
    hex_param_1 = (0.30557132498416676354695417221753673, 0.3466380372845119285587105)
    hex_param_2 = (0.30557132498416676354695417221753673, 1.44242681477454929484179805241)
    cut_start = (0, 1./sqrt(2))
    cut_end = (1./sqrt(2), 1./sqrt(2))

    def __init__(self, parent, packing_canvas, dash) :
        super().__init__(parent)
        self.default_canvas = QPointF(*self.to_canvas(self.default_param))
        self.hex_1_canvas = QPointF(*self.to_canvas(self.hex_param_1))
        self.hex_2_canvas = QPointF(*self.to_canvas(self.hex_param_2))
        self.cut_start_canvas = QPointF(*self.to_canvas(self.cut_start))
        self.cut_end_canvas = QPointF(*self.to_canvas(self.cut_end))
        self._is_dragging = False
        self._packing_canvas = packing_canvas
        self._packing_canvas.param = self.param
        self._dashboard = dash
    
    @property
    def packing_canvas(self) :
        return self._packing_canvas
    
    @property
    def param(self) :
        if not hasattr(self, '_param') :
            self._param = self.default_param
        return self._param

    def _reset(self) :
        self._param = self.default_param
        self._param_canvas = self.to_canvas(self._param)
        self.update()
        self._dashboard.param = self.param
        self._packing_canvas._reset()
        self._packing_canvas.param = self.param
    
    @property
    def param_canvas(self) :
        if not hasattr(self, '_param_canvas') :
            self._param_canvas = self.to_canvas(self.param)
        return self._param_canvas

    @property
    def param_bound_points(self) :
        if not hasattr(self, '_param_bound_points') :
            t = np.arange(0., self.max_t + 0.1, 0.1)
            s = (2*t)/(1 + 2*t**2)
            self._param_bound_points = np.array(list(map(self.to_canvas,zip(s,t))))
        return self._param_bound_points
 
    @property
    def unit_unif_points(self) :
        if not hasattr(self, '_unit_unif_points') :
            t = np.arange(0.00000000000001, self.max_t + 0.1, 0.1)
            s = (1 + 2*t**2 - sqrt(1 + 4*t**4))/(2*t)
            self._unit_unif_points = np.array(list(map(self.to_canvas,zip(s,t))))
        return self._unit_unif_points

    def _get_canvas_path(self, points) :
        (P1, P2) = bezier_control_points(points)
        path = QPainterPath()
        path.moveTo(QPointF(*points[0]))
        for p1,p2,p3 in zip(P1, P2, points[1:]) :
            path.cubicTo(*starmap(QPointF,(p1,p2,p3)))
        path.lineTo(QPointF(self.shift[0],points[-1][1]))
        return path 

    @property
    def param_bound_path(self) :
        if not hasattr(self, '_param_bound_path') :
            self._param_bound_path = self._get_canvas_path(self.param_bound_points)
        return self._param_bound_path
 
    @property
    def unit_unif_path(self) :
        if not hasattr(self, '_unit_unif_path') :
            self._unit_unif_path = self._get_canvas_path(self.unit_unif_points)
        return self._unit_unif_path

    @property
    def rect_unif_path(self) :
        if not hasattr(self, '_rect_unif_path') :
            self._rect_unif_path = QPainterPath()
            self._rect_unif_path.moveTo(self.cut_start_canvas)
            self._rect_unif_path.lineTo(self.cut_end_canvas)
        return self._rect_unif_path

    def draw_graphs(self, painter) :
        painter.setBrush(Qt.darkCyan)
        painter.fillPath(self.param_bound_path, painter.brush())
        painter.drawPath(self.unit_unif_path)
        painter.drawPath(self.rect_unif_path)

    def draw_params(self, painter) :
        painter.setBrush(Qt.black)
        painter.drawEllipse(self.default_canvas,4.,4.)
        painter.drawEllipse(self.hex_1_canvas,4.,4.)
        painter.drawEllipse(self.hex_2_canvas,4.,4.)
        painter.setBrush(Qt.red)
        painter.drawEllipse(QPointF(*self.param_canvas),4.,4.)
    
    def paintEvent(self, event) :
        painter = QPainter(self)
        painter.setRenderHint(QPainter.Antialiasing)
        self.draw_graphs(painter)
        self.draw_params(painter)

    def to_plane(self, p) :
        return ((p[0] - self.shift[0])/self.scale[0],
                (self.shift[1] - p[1])/self.scale[1])

    def to_canvas(self, p) :
        return (self.shift[0] + self.scale[0]*p[0],
                self.shift[1] - self.scale[1]*p[1])

    def mousePressEvent(self, event):
        if event.button() == Qt.LeftButton and not self._is_dragging:
            point = _get_point(event)
            diff = tuple(self.param_canvas[i] - point[i] for i in range(2))
            dist = sqrt(sum(d**2 for d in diff))
            if dist < self.drag_dist :
                self._is_dragging = True
                self._drag_diff = diff

    def _move_param(self, p, point = None) :
        if 0 < p[0] and p[1] <= self.max_t and p[0] <= (2*p[1])/(1 + 2*p[1]**2) :
            if point is not None :
                self._param_canvas = point
            else :
                self._param_canvas = self.to_canvas(p)
            self._param = p
            self.update()
            self._dashboard.param = self.param
            self._packing_canvas.param = self.param

    def mouseMoveEvent(self, event):
        if self._is_dragging :
            e_pt = _get_point(event)
            canv_pt = (e_pt[0]+self._drag_diff[0], e_pt[1]+self._drag_diff[1])
            p = self.to_plane(canv_pt)
            self._move_param(p, canv_pt)

    def mouseReleaseEvent(self, event):
        if event.button() == Qt.LeftButton and self._is_dragging :
            self._is_dragging = False

class ParamField(QLineEdit) :
    def __init__(self, dash) :
        super().__init__()
        self._dash = dash
        self.setFocusPolicy(Qt.ClickFocus)
        self.returnPressed.connect(self.process_on_return)

    def focusInEvent(self, event) :
        self.setReadOnly(False)
        super().focusInEvent(event)
        self._dash.paused = True

    def focusOutEvent(self, event) :
        super().focusInEvent(event)
        self.setReadOnly(True)
        self._dash.paused = False

    def process_on_return(self) :
        try :
            val = eval(self.text().replace('^','**'))
            if isinstance(val, float) :
                self._dash.field_update(val, self)
        except :
            self.clearFocus()

class Dashboard(QWidget) :
    
    tolerance = 0.0000001
    default_param = (sqrt(2.) - 1., 1/sqrt(2))

    def __init__(self, parent, packing) :
        super().__init__(parent)
        self._packing_canvas = packing
        # annoying width hack because of line elements
        self.setMaximumWidth(320)
        layout = QGridLayout(self)
        self._s_field = ParamField(self)
        self._t_field = ParamField(self)
        self._tr_one_label = QLabel()
        self._tr_two_label = QLabel()
        self._unif_label = QLabel()
        param_text = QLabel("Parameters (s,t) with 0 < s < 2t/(1 + 2t^2):")
        tr_one_text = QLabel("Trace of hol1 :")
        tr_two_text = QLabel("Trace of hol2 :")
        tr_unif_text = QLabel("Unif is ± arccosh(tr(hol2)/2)/arccosh(tr(hol1)/2) :")
        layout.addWidget(param_text,0,0,1,4)
        layout.addWidget(self._s_field,1,0,1,2)
        layout.addWidget(self._t_field,1,2,1,2)
        layout.addWidget(tr_one_text,2,0,1,4)
        layout.addWidget(self._tr_one_label,3,0,1,4)
        layout.addWidget(tr_two_text,4,0,1,4)
        layout.addWidget(self._tr_two_label,5,0,1,4)
        layout.addWidget(tr_unif_text,6,0,1,4)
        layout.addWidget(self._unif_label,7,0,1,4)
        self._reset_button = QPushButton("&Reset")
        self._reset_button.clicked.connect(self._reset)
        layout.addWidget(self._reset_button,8,0,1,1)
        self._dual_box = QCheckBox("&Duals")
        self._dual_box.stateChanged.connect(self._toggle_dual)
        layout.addWidget(self._dual_box,8,1,1,1)
        self._exp_box = QCheckBox("&Exp(cz)")
        self._exp_box.stateChanged.connect(self._toggle_exp)
        layout.addWidget(self._exp_box,8,2,1,1)
        self._color_box = QCheckBox("&Color")
        self._color_box.stateChanged.connect(self._toggle_color)
        layout.addWidget(self._color_box,8,3,1,1)
        # Color is on by default
        self._color_box.setCheckState(Qt.Checked)

    def _toggle_dual(self, state) :
        self._packing_canvas.show_dual = (state == Qt.Checked)
    
    def _toggle_exp(self, state) :
        self._packing_canvas.show_exp = (state == Qt.Checked)

    def _toggle_color(self, state) :
        self._packing_canvas.show_color = (state == Qt.Checked)

    def _reset(self) :
        # must be set after creation
        self._param_space._reset()

    def field_update(self, val, field) :
        if field is self._s_field :
            self._param_space._move_param((val, self.param[1]))
        elif field is self._t_field :
            self._param_space._move_param((self.param[0], val))

    @property
    def param(self) :
        if not hasattr(self, '_param') :
            self._param = self.default_param
        return self._param
            
    @param.setter
    def param(self, new) :
        self._s_field.clearFocus()
        self._t_field.clearFocus()
        self._param = new
        self.update()
    
    @property
    def paused(self) :
        if not hasattr(self,'_paused') :
            self._paused = False
        return self._paused
            
    @paused.setter
    def paused(self, new) :
        if isinstance(new, bool) :
            self._paused = new

    def paintEvent(self, event) :
        if not self.paused :
            self._tr_one = trace_one(*self.param)
            self._tr_two = trace_two(*self.param)
            if self.param != self.default_param :
                self._unif = tau_param(*self.param)
            else :
                self._unif = 1j
            
            self._s_field.setText('{:.10f}'.format(self.param[0]))
            self._t_field.setText('{:.10f}'.format(self.param[1]))
            self._tr_one_label.setText('{:.10f}'.format(self._tr_one))
            self._tr_two_label.setText('{:.10f}'.format(self._tr_two))
            self._unif_label.setText('{:.10f}'.format(self._unif))


from functools import partial

class PackingCanvas(QWidget) :
    
    default_param = (sqrt(2.) - 1., 1/sqrt(2))

    def __init__(self,parent) :
        super().__init__(parent)
        self._is_dragging  = False
        self._center_offset = (0.,0.)
        self._num_H = 5
        self._num_V  = 4
        self._zoom_buttons = [QPushButton("-",self), QPushButton("+",self)]
        self._H_buttons = [QPushButton("H -",self), QPushButton("H +",self)]
        self._V_buttons = [QPushButton("V -",self), QPushButton("V +",self)]
        for b in self._zoom_buttons + self._H_buttons + self._V_buttons :
            b.setFocusPolicy(Qt.NoFocus)
            b.setMaximumSize(30,30)
        for b in self._zoom_buttons :
            b.clicked.connect(self._zoom)
        for b in self._H_buttons :
            b.clicked.connect(self._H_delta)
        for b in self._V_buttons :
            b.clicked.connect(self._V_delta)

    def _reset(self) :
        self._center_offset = (0.,0.)
        del self._scale
        self._num_H = 5
        self._num_V = 4

    def _H_delta(self) :
        i = self._H_buttons.index(self.sender())
        new = self._num_H + 2*i - 1
        self._num_H = new if new >= 0 else 5
        # hack!
        self.param = self.param

    # gah, so repetative
    def _V_delta(self) :
        i = self._V_buttons.index(self.sender())
        new = self._num_V + 2*i - 1
        self._num_V = new if new >= 0 else 4
        # hack!
        self.param = self.param

    def _zoom(self) :
        i = self._zoom_buttons.index(self.sender())
        self.scale += 4.*i - 2.

    def resizeEvent(self, event) :
        super().resizeEvent(event)
        self.align_buttons()
    
    def align_buttons(self) :
        w,h  = self.width(), self.height()
        for i,b in enumerate(self._zoom_buttons) :
            b.move(w - 30., h - (i+1)*30.)
        for i,b in enumerate(self._H_buttons) :
            b.move(w - (i+1)*30., 0.)
        for i,b in enumerate(self._V_buttons) :
            b.move(w - (i+1)*30., 30.)

    @property
    def scale(self) :
        if not hasattr(self,'_scale') :
            self._scale = 30.
        return self._scale
            
    @scale.setter
    def scale(self,new) :
        self._scale = new if new > 0. else 30.
        self.update()

    @property
    def param(self) :
        if not hasattr(self,'_param') :
            self._param =  self.default_param
        return self._param
    
    @param.setter
    def param(self, new) :
        self._param = new
        s_right = shear_right(*self.param)
        s_up = shear_up(*self.param)
        self._data = fern(self._num_H,self._num_V,s_right,s_up)
        # branching will be a problem
        x = None
        if self.show_exp and self.param != self.default_param :
            x = sqrt(conj_param_squared(*self.param))
        self._circles = all_circles(self._data, self.show_dual, x)
        self.update()

    def _update_from_toggle(self) :
        x = None
        if self.show_exp and self.param != self.default_param :
            x = sqrt(conj_param_squared(*self.param))
        self._circles = all_circles(self._data, self.show_dual, x)
        self.update()

    @property
    def show_dual(self) :
        if not hasattr(self,'_show_dual') :
            self._show_dual = False
        return self._show_dual

    @show_dual.setter
    def show_dual(self, new) :
        if self.show_dual != new :
            self._show_dual = new
            self._update_from_toggle()

    @property
    def show_exp(self) :
        if not hasattr(self,'_show_exp') :
            self._show_exp = False
        return self._show_exp
    
    @show_exp.setter
    def show_exp(self, new) :
        if self.show_exp != new :
            self._show_exp = new
            self._update_from_toggle()

    @property
    def show_color(self) :
        if not hasattr(self,'_show_color') :
            self._show_color = False
        return self._show_color
    
    @show_color.setter
    def show_color(self, new) :
        if self.show_color != new :
            self._show_color = new
            self.update()

    @property
    def circles(self) :
        return self._circles
    
    def paintEvent(self, event) :
        super().paintEvent(event)
        painter = QPainter(self)
        painter.setRenderHint(QPainter.Antialiasing)
        self.draw_circles(painter)

    def draw_circles(self, painter) :
        
        c_x = self._center_offset[0] + self.width()/2.
        c_y = self._center_offset[1] + self.height()/2.
        
        for type, circles in self.circles.items() :
            
            if self.show_color :
                if type == 'small' :
                    painter.setBrush(QColor(255, 0, 0, 20))
                elif type == 'large' :
                    painter.setBrush(QColor(0, 0, 255, 20))
                else :
                    painter.setBrush(QColor(0, 0, 0, 0))
            
            for c in circles :
                if isinstance(c[1],complex) :
                    z1 = self.scale*c[0]
                    z2 = self.scale*c[1]
                    v = (z2-z1)/abs(z2-z1)
                    s = 2*max(c_x,c_y)
                    p1 = z1 + s*v
                    p2 = z1 - s*v
                    l = QLineF(c_x + p1.real, c_y - p1.imag,
                               c_x + p2.real, c_y - p2.imag)
                    painter.drawLine(l)
                else :
                    r = self.scale*c[1] # scale factor
                    # if r > 10000. : continue
                    x = self.scale*c[0].real
                    y = self.scale*c[0].imag
                    painter.drawEllipse(QPointF(c_x+x,c_y-y), r, r)

    def mousePressEvent(self, event) :
        if event.button() == Qt.LeftButton and not self._is_dragging:
            self._is_dragging = True
            self._drag_start = _get_point(event)
            self._ref_offset = self._center_offset

    def mouseMoveEvent(self, event) :
        if self._is_dragging :
            p = _get_point(event)
            new_x = self._ref_offset[0] + p[0] - self._drag_start[0]
            new_y = self._ref_offset[1] + p[1] - self._drag_start[1]
            self._center_offset = (new_x,new_y)
            self.update()

    def mouseReleaseEvent(self, event):
        if event.button() == Qt.LeftButton and self._is_dragging :
            self._is_dragging = False

class ApplicationWindow(QMainWindow) :
    
    def __init__(self):
        super().__init__()
        self.setAttribute(Qt.WA_DeleteOnClose)
        self.setWindowTitle("Circle Packing")
    
        central = QWidget(self)
        self.setCentralWidget(central)
        layout = QGridLayout(central)
        self._packing = PackingCanvas(central)
        self._dash = Dashboard(central, self._packing)
        self._param = ParamSpace(central, self._packing, self._dash)
        self._dash._param_space = self._param
        layout.addWidget(self._packing, 0, 0, 3, 3)
        layout.addWidget(self._dash,2,3,1,1)
        layout.addWidget(self._param, 0, 3, 2, 1)
    
    def fileQuit(self):
        self.close()
    
    def closeEvent(self, ce):
        self.fileQuit()
    
    def about(self):
        QMessageBox.about(self, "About", "Stuff")

if __name__ == '__main__':
    import sys
    app = QApplication(sys.argv)
    app_window = ApplicationWindow()
    
    app_window.showMaximized()
    sys.exit(app.exec_())
