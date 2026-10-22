#!/usr/bin/env python3
"""Compact geometry and exhaustive local-support diagrams for the 3D notes.

The data come from rvb.octahedral; these are geometric/support diagrams,
not a numerical correlation calculation or a statement about singlet phases.
"""
from __future__ import annotations

from itertools import combinations, product
from pathlib import Path
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import FancyArrowPatch, Polygon
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from rvb.octahedral import LOCAL_STATES, binary_rank

OUT = ROOT / "results" / "collaborator-notes"
GREEN = "#22866f"
PURPLE = "#8b4b9b"
BLUE = "#2675a6"
ORANGE = "#d78c29"
GRAY = "#c8ced3"
INK = "#253746"
PORTS = np.array(((1,0,0),(-1,0,0),(0,1,0),(0,-1,0),(0,0,1),(0,0,-1)),float)
# An oblique, nondegenerate projection shared by every octahedron.
PROJECT = np.array(((.88,-.59,.08),(-.31,-.46,.92)))
EDGES = [(a,b) for a,b in combinations(range(6),2) if a//2 != b//2]


def project(points):
    return np.asarray(points) @ PROJECT.T


def wire(ax, points, color=GRAY, lw=.55, alpha=1, zorder=1):
    xy = project(points)
    for a,b in EDGES:
        ax.plot(*xy[[a,b]].T,color=color,lw=lw,alpha=alpha,zorder=zorder)
    return xy


def save(fig,name):
    OUT.mkdir(parents=True,exist_ok=True)
    fig.savefig(OUT/f"{name}.pdf",facecolor="white")
    fig.savefig(OUT/f"{name}.png",dpi=220,facecolor="white")
    plt.close(fig)


def catalogue():
    fig,axs=plt.subplots(4,4,figsize=(6.7,3.65))
    fig.subplots_adjust(left=.015,right=.985,bottom=.09,top=.99,hspace=.13,wspace=.02)
    axis_order=((),(0,),(1,),(2,),(0,1),(0,2),(1,2),(0,1,2))
    ordered=[next(s for s in LOCAL_STATES if not s.is_trimer and s.dimer_axes==a)
             for a in axis_order]
    ordered += [next(s for s in LOCAL_STATES if s.face_signs==signs)
                for signs in product((1,-1),repeat=3)]
    supports={s.support for s in ordered}
    assert len(supports)==16
    # Six local bits, two independent homogeneous equations: dimension four.
    # This local vector space is distinct from the globally constrained
    # half-edge ownership fields, which form an affine space because every
    # shared site has complementary owners at its two ends.
    parity_checks=(0b001111,0b110011)
    assert binary_rank(parity_checks)==2
    assert all((a ^ b) in supports for a in supports for b in supports)
    assert supports=={a for a in range(64)
                     if all((a & row).bit_count()%2==0 for row in parity_checks)}
    for ax,state in zip(axs.flat,ordered):
        incoming=np.array(state.incoming,dtype=bool)
        parity=[int(incoming[2*i]) ^ int(incoming[2*i+1]) for i in range(3)]
        assert parity[0]==parity[1]==parity[2]
        xy=wire(ax,PORTS)
        color=PURPLE if state.is_trimer else GREEN
        if state.is_trimer:
            face=xy[incoming]
            ax.add_patch(Polygon(face,facecolor=PURPLE,alpha=.15,edgecolor="none",zorder=2))
            ax.plot(*np.vstack((face,face[0])).T,color=color,lw=1.35,zorder=3)
            signs="".join("+" if s>0 else "-" for s in state.face_signs)
            label=r"$T_{"+signs+r"}$"
        else:
            for axis in state.dimer_axes:
                ax.plot(*xy[2*axis:2*axis+2].T,color=color,lw=1.55,zorder=3)
            label="Empty" if not state.dimer_axes else "$"+"+".join("D_{"+"xyz"[a]+"}" for a in state.dimer_axes)+"$"
        ax.scatter(*xy.T,c=np.where(incoming,color,"#e3e6e9"),s=15,
                   edgecolors=np.where(incoming,color,"#aeb7be"),linewidths=.5,zorder=4)
        # External radial spokes keep the six ownership arrows separate from
        # all physical dimer chords, face edges and the octahedral wireframe.
        for p,point in enumerate(xy):
            unit=point/np.linalg.norm(point)
            near=point+.12*unit
            far=point+.59*unit
            start,end=(far,near) if incoming[p] else (near,far)
            ax.add_patch(FancyArrowPatch(start,end,arrowstyle="-|>",
                                        mutation_scale=6.4,lw=.85,
                                        color=color if incoming[p] else "#687781",
                                        shrinkA=0,shrinkB=0,zorder=5))
        ax.text(0,1.79,label,fontsize=8.2,ha="center",va="center",color=INK)
        if not state.dimer_axes and not state.is_trimer:
            for label,position in [("+x",(1.60,-.30)),("+y",(-1.29,-1.09)),("+z",(-.35,1.14))]:
                ax.text(*position,label,fontsize=5.8,
                        color="#596873",ha="center",va="bottom")
        ax.set(xlim=(-2.05,2.05),ylim=(-1.62,2.02),aspect="equal")
        ax.axis("off")
    fig.text(.5,.043,"Colored arrows point inward: owned here.  Gray arrows point outward: owned by a neighbour.",
             fontsize=7.1,ha="center",color=INK)
    fig.text(.5,.012,"Arrows show ownership on all six cubic ports.  Trimer signs are ordered (x, y, z).",
             fontsize=7.1,ha="center",color=INK)
    save(fig,"octahedral-catalogue")


def geometry():
    fig,axs=plt.subplots(1,2,figsize=(6.7,2.12),gridspec_kw={"width_ratios":[1,1.15]})
    fig.subplots_adjust(left=.015,right=.985,bottom=.10,top=.89,wspace=.08)
    ax=axs[0]
    tetra=np.array(((1,1,1),(1,-1,-1),(-1,1,-1),(-1,-1,1)),float)
    txy=project(tetra)
    for a,b in combinations(range(4),2):
        ax.plot(*txy[[a,b]].T,color="#8396a2",lw=.7,ls=(0,(3,2)),zorder=1)
    # One old-face triangle, and one old-vertex-star triangle. There are four of each.
    old_face=np.array(((tetra[0]+tetra[1])/2,(tetra[0]+tetra[2])/2,(tetra[1]+tetra[2])/2))
    old_star=np.array(((tetra[0]+tetra[1])/2,(tetra[0]+tetra[2])/2,(tetra[0]+tetra[3])/2))
    ax.add_patch(Polygon(project(old_face),fc=BLUE,alpha=.22,ec="none",zorder=2))
    ax.add_patch(Polygon(project(old_star),fc=ORANGE,alpha=.25,ec="none",zorder=2))
    xy=wire(ax,PORTS,color=BLUE,lw=.7,zorder=3)
    ax.scatter(*xy.T,s=15,c=BLUE,zorder=4)
    ax.scatter(*txy.T,s=8,c="#8396a2",zorder=4)
    ax.text(.02,1.04,"Tetrahedron → inscribed octahedron",transform=ax.transAxes,
            ha="left",fontsize=9,color=INK)
    ax.text(1.08,.65,"4 old-face\ntriangles",fontsize=7,color=BLUE)
    ax.text(1.08,-.12,"4 vertex-star\ntriangles",fontsize=7,color="#a36515")
    ax.text(.02,-.07,"Six edge midpoints become the six covering sites.",
            transform=ax.transAxes,fontsize=7,color=INK)
    ax.set(xlim=(-1.65,3.2),ylim=(-1.6,1.6),aspect="equal")
    ax.axis("off")

    ax=axs[1]
    centres=np.array(list(product((0,1),repeat=3)),float)
    shared=[]
    for c in centres:
        wire(ax,c+.5*PORTS,color="#c7ced4",lw=.45,alpha=.8,zorder=1)
        for d in np.eye(3):
            if np.all(c+d<=1):
                edge=project([c,c+d])
                ax.plot(*edge.T,color=ORANGE,lw=1.1,zorder=2)
                shared.append(c+.5*d)
    cxy=project(centres)
    mxy=project(shared)
    ax.scatter(*cxy.T,c=ORANGE,s=13,zorder=4)
    ax.scatter(*mxy.T,c=BLUE,s=15,zorder=5)
    ax.text(.01,1.04,"Corner-sharing octahedra; cubic centres",transform=ax.transAxes,
            ha="left",fontsize=9,color=INK)
    ax.text(.01,-.07,"Every blue midpoint has exactly two owners.",
            transform=ax.transAxes,fontsize=7,color=INK)
    ax.legend(handles=[Line2D([],[],color=ORANGE,marker="o",ms=3,lw=.8,label="centres / bonds"),
                       Line2D([],[],color=BLUE,marker="o",ms=3,lw=0,label="physical sites")],
              loc="lower right",bbox_to_anchor=(1.02,-.03),fontsize=6.3,
              frameon=False,handlelength=1.2,labelspacing=.3,borderpad=0)
    ax.set(xlim=(-1.20,2.0),ylim=(-1.37,1.64),aspect="equal")
    ax.axis("off")
    save(fig,"geometry-origin")


if __name__=="__main__":
    plt.rcParams.update({"font.family":"DejaVu Sans","pdf.fonttype":42,"ps.fonttype":42})
    catalogue()
    geometry()
