# -*- coding: utf-8 -*-
"""
svgplot.py —— 极简 SVG 绘图库（纯标准库，无 matplotlib / numpy 依赖）
用途：出期刊级矢量图（SVG 可直接投稿，或用 Inkscape / 浏览器转 PDF）
"""
import math

FONT = "Arial, Helvetica, sans-serif"

# Okabe-Ito 色盲友好配色
C_BLUE   = "#0072B2"
C_ORANGE = "#E69F00"
C_GREEN  = "#009E73"
C_RED    = "#D55E00"
C_PURPLE = "#CC79A7"
C_SKY    = "#56B4E9"
C_YELLOW = "#F0E442"
C_GRAY   = "#7F7F7F"
C_DARK   = "#222222"


def esc(s):
    return (str(s).replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;"))


class Fig:
    def __init__(self, w, h, bg="#ffffff"):
        self.w, self.h, self.bg = w, h, bg
        self.o = []

    # ---------- 基本图元 ----------
    def text(self, x, y, s, size=7, anchor="middle", fill=C_DARK, weight="normal",
             style="normal", rotate=None, baseline="middle", opacity=None, halo=False):
        dy = 0.0
        if baseline == "middle":
            dy = size * 0.35
        elif baseline == "top":
            dy = size * 0.90
        tr = ' transform="rotate(%.2f %.2f %.2f)"' % (rotate, x, y) if rotate else ""
        op = ' fill-opacity="%.2f"' % opacity if opacity is not None else ""
        hl = (' stroke="#ffffff" stroke-width="2.0" stroke-linejoin="round"'
              ' paint-order="stroke"') if halo else ""
        self.o.append(
            '<text x="%.2f" y="%.2f" font-family="%s" font-size="%.2f" fill="%s" '
            'font-weight="%s" font-style="%s" text-anchor="%s"%s%s%s>%s</text>'
            % (x, y + dy, FONT, size, fill, weight, style, anchor, tr, op, hl, esc(s)))

    def line(self, x1, y1, x2, y2, stroke="#333333", w=1.0, dash=None, cap="butt",
             opacity=None):
        d = ' stroke-dasharray="%s"' % dash if dash else ""
        op = ' stroke-opacity="%.2f"' % opacity if opacity is not None else ""
        self.o.append(
            '<line x1="%.2f" y1="%.2f" x2="%.2f" y2="%.2f" stroke="%s" '
            'stroke-width="%.2f" stroke-linecap="%s"%s%s/>'
            % (x1, y1, x2, y2, stroke, w, cap, d, op))

    def rect(self, x, y, w, h, fill="none", stroke="none", sw=0.9, dash=None,
             rx=0, opacity=None, sopacity=None):
        d = ' stroke-dasharray="%s"' % dash if dash else ""
        op = ' fill-opacity="%.2f"' % opacity if opacity is not None else ""
        so = ' stroke-opacity="%.2f"' % sopacity if sopacity is not None else ""
        self.o.append(
            '<rect x="%.2f" y="%.2f" width="%.2f" height="%.2f" rx="%.2f" '
            'fill="%s" stroke="%s" stroke-width="%.2f"%s%s%s/>'
            % (x, y, w, h, rx, fill, stroke, sw, d, op, so))

    def circle(self, cx, cy, r, fill="none", stroke="none", sw=1.0, opacity=None):
        op = ' fill-opacity="%.2f"' % opacity if opacity is not None else ""
        self.o.append(
            '<circle cx="%.2f" cy="%.2f" r="%.2f" fill="%s" stroke="%s" '
            'stroke-width="%.2f"%s/>' % (cx, cy, r, fill, stroke, sw, op))

    def poly(self, pts, fill="none", stroke="none", sw=1.0, opacity=None):
        p = " ".join("%.2f,%.2f" % (a, b) for a, b in pts)
        op = ' fill-opacity="%.2f"' % opacity if opacity is not None else ""
        self.o.append(
            '<polygon points="%s" fill="%s" stroke="%s" stroke-width="%.2f"%s/>'
            % (p, fill, stroke, sw, op))

    def polyline(self, pts, stroke="#333333", sw=1.0, fill="none", dash=None):
        p = " ".join("%.2f,%.2f" % (a, b) for a, b in pts)
        d = ' stroke-dasharray="%s"' % dash if dash else ""
        self.o.append(
            '<polyline points="%s" fill="%s" stroke="%s" stroke-width="%.2f" '
            'stroke-linejoin="round" stroke-linecap="round"%s/>' % (p, fill, stroke, sw, d))

    def path(self, d, fill="none", stroke="none", sw=1.0, dash=None):
        ds = ' stroke-dasharray="%s"' % dash if dash else ""
        self.o.append('<path d="%s" fill="%s" stroke="%s" stroke-width="%.2f"%s/>'
                      % (d, fill, stroke, sw, ds))

    def arrow(self, x1, y1, x2, y2, stroke=C_GRAY, w=0.9, dash=None, head=3.2):
        """直线箭头（末端三角）"""
        ang = math.atan2(y2 - y1, x2 - x1)
        bx, by = x2 - head * math.cos(ang), y2 - head * math.sin(ang)
        self.line(x1, y1, bx, by, stroke=stroke, w=w, dash=dash)
        p1 = (x2, y2)
        p2 = (bx - head * 0.45 * math.sin(ang), by + head * 0.45 * math.cos(ang))
        p3 = (bx + head * 0.45 * math.sin(ang), by - head * 0.45 * math.cos(ang))
        self.poly([p1, p2, p3], fill=stroke)

    # ---------- 输出 ----------
    def render(self):
        return ('<svg xmlns="http://www.w3.org/2000/svg" width="%.1fpt" height="%.1fpt" '
                'viewBox="0 0 %.1f %.1f" version="1.1">\n'
                '<rect x="0" y="0" width="%.1f" height="%.1f" fill="%s"/>\n%s\n</svg>'
                % (self.w, self.h, self.w, self.h, self.w, self.h, self.bg, "\n".join(self.o)))

    def save(self, path):
        with open(path, "w", encoding="utf-8") as fh:
            fh.write('<?xml version="1.0" encoding="UTF-8"?>\n' + self.render())


class Ax:
    """直角坐标系：box=(left, top, right, bottom)，y 轴向下"""
    def __init__(self, fig, box, xlim, ylim, xlog=True, ylog=True):
        self.f = fig
        self.x0, self.y0, self.x1, self.y1 = box
        self.xlim, self.ylim = xlim, ylim
        self.xlog, self.ylog = xlog, ylog

    def _t(self, v, log):
        return math.log10(v) if log else v

    def X(self, v):
        a, b = self._t(self.xlim[0], self.xlog), self._t(self.xlim[1], self.xlog)
        return self.x0 + (self._t(v, self.xlog) - a) / (b - a) * (self.x1 - self.x0)

    def Y(self, v):
        a, b = self._t(self.ylim[0], self.ylog), self._t(self.ylim[1], self.ylog)
        return self.y1 - (self._t(v, self.ylog) - a) / (b - a) * (self.y1 - self.y0)

    def frame(self, left=True, bottom=True, top=False, right=False, color="#333333", w=0.9):
        f = self.f
        if left:
            f.line(self.x0, self.y0, self.x0, self.y1, stroke=color, w=w)
        if bottom:
            f.line(self.x0, self.y1, self.x1, self.y1, stroke=color, w=w)
        if top:
            f.line(self.x0, self.y0, self.x1, self.y0, stroke=color, w=w)
        if right:
            f.line(self.x1, self.y0, self.x1, self.y1, stroke=color, w=w)

    def grid(self, xs=None, ys=None, color="#EBEBEB", w=0.6):
        f = self.f
        for v in (xs or []):
            f.line(self.X(v), self.y0, self.X(v), self.y1, stroke=color, w=w)
        for v in (ys or []):
            f.line(self.x0, self.Y(v), self.x1, self.Y(v), stroke=color, w=w)

    def xticks(self, vals, labels=None, size=6.5, ticklen=2.5, pad=6,
               color="#333333", weight="normal", rotate=None):
        f = self.f
        for i, v in enumerate(vals):
            x = self.X(v)
            f.line(x, self.y1, x, self.y1 + ticklen, stroke=color, w=0.8)
            lab = labels[i] if labels else ("%g" % v)
            if lab == "":
                continue
            if rotate:
                f.text(x, self.y1 + pad + 1.5, lab, size=size, anchor="end",
                       rotate=rotate, fill="#333333", weight=weight)
            else:
                f.text(x, self.y1 + pad + 1.5, lab, size=size, anchor="middle",
                       fill="#333333", weight=weight)

    def yticks(self, vals, labels=None, size=6.5, ticklen=2.5, pad=5, color="#333333"):
        f = self.f
        for i, v in enumerate(vals):
            y = self.Y(v)
            f.line(self.x0 - ticklen, y, self.x0, y, stroke=color, w=0.8)
            lab = labels[i] if labels else ("%g" % v)
            if lab == "":
                continue
            f.text(self.x0 - pad, y, lab, size=size, anchor="end", fill="#333333")

    def xlabel(self, s, size=7.5, dy=17):
        self.f.text((self.x0 + self.x1) / 2, self.y1 + dy, s, size=size)

    def ylabel(self, s, size=7.5, dx=16):
        self.f.text(self.x0 - dx, (self.y0 + self.y1) / 2, s, size=size, rotate=-90)
