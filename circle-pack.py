#!/usr/bin/env python3

from PyQt5.QtCore import (Qt, QPointF, QLineF, QSize)
from PyQt5.QtGui import (QPainter, QPainterPath, QValidator, QColor)
from PyQt5.QtWidgets import (QApplication, QMainWindow, QWidget,
                             QPushButton, QCheckBox, QLabel,
                             QLineEdit, QGridLayout)

from itertools import starmap
from circletools import *

class ParamSpace(QWidget) :

    err_dist = 6
    max_t = 2.5
    scale = (400.,195.)
    shift = (20.,490.)
    default_param = (sqrt(2.) - 1., 1./sqrt(2))
    hex_params = [(0.30557132498416676354695417221753673, 0.3466380372845119285587105),
                  (0.30557132498416676354695417221753673, 1.44242681477454929484179805241)]
    fun_params = [(0.412,0.6485),(0.40088,0.5691), (0.580306,1./sqrt(2)), (0.45616,0.58936), (0.25,1.0), (0.40054,1./sqrt(2))]
    cut_start = (0, 1./sqrt(2))
    cut_end = (1./sqrt(2), 1./sqrt(2))

    def __init__(self, parent, packing_canvas, dash) :
        super().__init__(parent)
        self.default_canvas = self.to_canvas(self.default_param)
        self.hex_canvas = [self.to_canvas(x) for x in self.hex_params]
        self.fun_canvas = [self.to_canvas(x) for x in self.fun_params]
        self.special_canvas = self.fun_canvas + self.hex_canvas
        self.cut_start_canvas = self.to_canvas(self.cut_start)
        self.cut_end_canvas = self.to_canvas(self.cut_end)
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
            self._param_bound_points = np.array(list(map(lambda p : (p.x(),p.y()),
                                                     map(self.to_canvas,zip(s,t)))))
        return self._param_bound_points
 
    @property
    def unit_unif_points(self) :
        if not hasattr(self, '_unit_unif_points') :
            t = np.arange(0.00000000000001, self.max_t + 0.1, 0.1)
            s = (1 + 2*t**2 - sqrt(1 + 4*t**4))/(2*t)
            self._unit_unif_points = np.array(list(map(lambda p : (p.x(),p.y()),
                                                   map(self.to_canvas,zip(s,t)))))
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
        for p in self.hex_canvas :
            painter.drawEllipse(p,4.,4.)
        painter.setBrush(Qt.blue)
        for p in self.fun_canvas :
            painter.drawEllipse(p,4.,4.) 
        painter.setBrush(Qt.red)
        painter.drawEllipse(self.param_canvas,4.,4.)
    
    def paintEvent(self, event) :
        painter = QPainter(self)
        painter.setRenderHint(QPainter.Antialiasing)
        self.draw_graphs(painter)
        self.draw_params(painter)

    def to_plane(self, p) :
        return ((p.x() - self.shift[0])/self.scale[0],
                (self.shift[1] - p.y())/self.scale[1])

    def to_canvas(self, p) :
        return QPointF(self.shift[0] + self.scale[0]*p[0],
                       self.shift[1] - self.scale[1]*p[1])

    def mousePressEvent(self, event):
        if event.button() == Qt.LeftButton and not self._is_dragging:
            point = event.pos()
            diff = self.param_canvas - point
            dist = QPointF.dotProduct(diff,diff) 
            if dist < self.err_dist :
                self._is_dragging = True
                self._drag_diff = diff
            else :
                self._jump_click = True;

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

    def snap_special_or_self(self, p) :
        for q in self.special_canvas :
            diff = q - p
            dist = QPointF.dotProduct(diff,diff)
            if dist < self.err_dist :
                return q
        return p

    def mouseMoveEvent(self, event):
        if self._is_dragging :
            e_pt = event.pos()
            canv_pt = e_pt + self._drag_diff
            p = self.to_plane(canv_pt)
            self._move_param(p, canv_pt)

    def mouseReleaseEvent(self, event):
        if event.button() == Qt.LeftButton :
            if self._is_dragging :
                self._is_dragging = False
            elif self._jump_click :
                e_pt = event.pos()
                canv_pt = self.snap_special_or_self(e_pt)
                p = self.to_plane(canv_pt)
                self._move_param(p, canv_pt)
                self._jump_click = False

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
        self._jump_click  = False
        self._center_offset = QPointF(0.,0.)
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
        self._center_offset = QPointF(0.,0.)
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
        
        c_x = self._center_offset.x() + self.width()/2.
        c_y = self._center_offset.y() + self.height()/2.
        
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
            self._drag_start = event.pos()
            self._ref_offset = self._center_offset

    def mouseMoveEvent(self, event) :
        if self._is_dragging :
            p = event.pos()
            self._center_offset = self._ref_offset + p - self._drag_start
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
