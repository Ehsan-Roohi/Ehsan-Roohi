"""Rebuild the English lecture from the method guide and measured plots.
Requires reportlab and pypdf. Does not rerun numerical simulations.
"""
from pathlib import Path
import hashlib, html, json, re
from reportlab.pdfgen import canvas
from reportlab.platypus import SimpleDocTemplate, Paragraph, Spacer, Table, TableStyle, PageBreak, Image, KeepTogether, Flowable
from reportlab.lib import colors
from reportlab.lib.styles import ParagraphStyle
from reportlab.lib.pagesizes import A4
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
from pypdf import PdfReader

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'output/pdf';OUT.mkdir(parents=True,exist_ok=True)
PDF=OUT/'CFD_Algorithms_Lecture_2026.pdf'
BLUE=colors.HexColor('#12344D');TEAL=colors.HexColor('#087F8C');LIGHT=colors.HexColor('#EDF5F7');GREY=colors.HexColor('#52616B')
W,H=A4;WIDTH=W-88
for name,file in [('Lecture','arial.ttf'),('LectureBold','arialbd.ttf'),('LectureItalic','ariali.ttf')]:
    pdfmetrics.registerFont(TTFont(name,str(Path('C:/Windows/Fonts')/file)))
pdfmetrics.registerFontFamily('Lecture',normal='Lecture',bold='LectureBold',italic='LectureItalic',boldItalic='LectureBold')
S={
 'body':ParagraphStyle('body',fontName='Lecture',fontSize=10.2,leading=14.3,textColor=BLUE,spaceAfter=8),
 'small':ParagraphStyle('small',fontName='Lecture',fontSize=8.6,leading=12,textColor=GREY,spaceAfter=5),
 'title':ParagraphStyle('title',fontName='LectureBold',fontSize=24,leading=29,textColor=BLUE,spaceAfter=13),
 'sub':ParagraphStyle('sub',fontName='LectureBold',fontSize=12.5,leading=17,textColor=TEAL,spaceBefore=8,spaceAfter=7),
 'card':ParagraphStyle('card',fontName='LectureBold',fontSize=11,leading=15,textColor=BLUE,spaceAfter=6)}
story=[];plan=[];seen=[]
def p(text,style='body'):return Paragraph(text,S[style])
def add(text,style='body'):story.append(p(text,style))
def heading(kicker,title):
    if plan:story.append(PageBreak())
    plan.append(title);add(kicker.upper(),'small');add(title,'title')
def sub(title):add(title,'sub')
def bullets(items):
    for item in items:add('• '+item)
def box(label,text):
    t=Table([[p(label,'card')],[p(text)]],colWidths=[WIDTH])
    t.setStyle(TableStyle([('BACKGROUND',(0,0),(-1,-1),LIGHT),('BOX',(0,0),(-1,-1),.6,TEAL),('LEFTPADDING',(0,0),(-1,-1),10),('RIGHTPADDING',(0,0),(-1,-1),10),('TOPPADDING',(0,0),(-1,0),8),('BOTTOMPADDING',(0,-1),(-1,-1),5)]))
    story.extend([t,Spacer(1,10)])
def eq(text):box('Key equation',text)
def fig(name,caption,maxheight=400):
    image=Image(str(ROOT/'figures'/name));r=min(WIDTH/image.imageWidth,maxheight/image.imageHeight)
    image.drawWidth=image.imageWidth*r;image.drawHeight=image.imageHeight*r
    story.extend([image,Spacer(1,7),p(caption,'small')])
def table(rows,widths):
    t=Table([[p(str(v),'small') for v in row] for row in rows],colWidths=widths,repeatRows=1,hAlign='LEFT')
    t.setStyle(TableStyle([('BACKGROUND',(0,0),(-1,0),LIGHT),('VALIGN',(0,0),(-1,-1),'TOP'),('LINEBELOW',(0,0),(-1,0),1,TEAL),('LINEBELOW',(0,1),(-1,-1),.3,colors.HexColor('#D8E4E8')),('LEFTPADDING',(0,0),(-1,-1),6),('RIGHTPADDING',(0,0),(-1,-1),6),('TOPPADDING',(0,0),(-1,-1),7),('BOTTOMPADDING',(0,0),(-1,-1),7)]))
    story.extend([t,Spacer(1,10)])
def prose(text):
    d={r'\widehat F':'Fhat',r'\lambda':'λ',r'\alpha':'α',r'\beta':'β',r'\varepsilon':'ε',r'\tau':'τ',r'\rho':'ρ',r'\vert':'|',r'\pm':'±',r'\max':'max',r'\min':'min',r'\Delta':'Δ',r'\,':' '}
    for a,b in d.items():text=text.replace(a,b)
    text=text.replace('$','').replace('{','').replace('}','').replace('\\','')
    text=html.escape(text).replace('—','-').replace('–','-')
    return re.sub(r'\*\*(.*?)\*\*',r'<b>\1</b>',text)
methods={}
for line in (ROOT/'ALGORITHM_GUIDE.md').read_text(encoding='utf-8').splitlines():
    if line.startswith('| `'):
        cols=line.strip('|').split('|');assert len(cols)==3
        key=re.search(r'`([^`]+)`',cols[0])[1]
        methods[key]=tuple(prose(c.replace('`','').strip()) for c in cols)
ADV={
'godunov':'A direct benchmark for approximate fluxes; retains the complete local wave structure.',
'roe':'A compact wave-by-wave treatment with explicit contact and acoustic information.',
'roe-nc':'Isolates the effect of Roe entropy correction in a diagnostic experiment.',
'hlle':'Few intermediate-state calculations; useful as the implemented HLLC fallback.',
'hllc':'Restores the contact wave without an iterative exact star-pressure solve.',
'rusanov':'Short, transparent implementation; useful for debugging and comparison.',
'global-lf':'A single speed makes the chosen dissipation easy to inspect.',
'vanleer':'Directional physical fluxes with a smooth subsonic-to-supersonic construction.',
'steger-warming':'A systematic eigenvalue split provides a clear upwind construction.',
'ausm':'Separates mass transport from pressure, making family changes easy to trace.',
'ausm-plus':'Improves Mach/pressure splitting without a local iterative Riemann solve.',
'ausm-up':'Adds explicit pressure/velocity coupling for tests across Mach regimes.',
'ausm-up2':'Tests another pressure-dissipation mechanism with the same mass-correction backbone.',
'slau2':'A distinct all-speed splitting design for controlled comparisons.',
'ec-lf':'Separates an entropy-conservative central construction from added dissipation.',
'first':'A minimal conservative baseline for verifying higher-order modifications.',
'muscl':'A compact second-order stencil balancing centered slopes and jump control.',
'muscl-minmod':'Simple, cautious slopes that are easy to verify by hand.',
'muscl-vanleer':'A smooth harmonic transition between permitted slopes.',
'muscl-superbee':'Can keep monotone interfaces more compact than gentler limiters.',
'muscl-vanalbada':'A smooth rational response without abrupt one-sided slope selection.',
'eno2':'Avoids the larger local jump by choosing the less troubled stencil.',
'cweno3':'One cell polynomial supports faces and geometric-source integration.',
'cweno5':'Higher smooth accuracy while retaining a cell-based polynomial.',
'cweno3-char':'Weights individual waves instead of mixing conservative components.',
'cweno5-char':'Combines wave weighting with a higher-order smooth polynomial.',
'weno3-js':'A small characteristic face stencil with nonlinear smoothness weighting.',
'weno3-z':'Improves recovery of optimal smooth weights relative to this JS construction.',
'weno5-js':'Fifth-order candidate blending when smooth stencil data are available.',
'weno5-z':'Uses a global smoothness difference to reduce unnecessary weight distortion.',
'weno7-js':'Raises formal face accuracy on a sufficiently smooth wider stencil.',
'teno3':'Kept stencils recover renormalized optimal weights instead of continuously altered weights.',
'teno5':'Targeted detection rejects troubled fifth-order candidate stencils.',
'teno7':'Tests targeted selection on a wider, higher-order face stencil.',
'muscl-thinc-bvd':'A bounded nonpolynomial candidate represents part of an interface inside a cell.'}
LIMIT={
'roe-nc':'Omitting the entropy correction can admit nonphysical expansion behavior; this diagnostic is not a general replacement for corrected Roe.',
'global-lf':'A fast region sets damping everywhere, increasing diffusion in slow regions. This is the semidiscrete flux, not a classical fully discrete LF update.',
'ausm':'The original split lacks the dedicated pressure/velocity corrections of AUSM+-up; all-Mach behavior needs separate tests.',
'ausm-up2':'Pressure dissipation and the reference-Mach setting can influence diffusion. The Sod study does not establish superiority across Mach regimes.',
'slau2':'Changing the dissipation design does not guarantee smaller errors in every flow. These shock-tube runs do not validate a low-Mach asymptotic limit.',
'ec-lf':'LF dissipation can spread contacts. A checked central face entropy identity is not a proof for the full reconstructed and time-integrated scheme.',
'muscl-vanleer':'Remains second order in smooth regions and flattens at extrema under the implemented same-sign condition. Contacts can still spread.',
'muscl-vanalbada':'Remains second order and returns zero at opposing slope signs; smooth extrema and transported contacts can lose resolution.',
'cweno3':'Candidate fitting and nonlinear weighting cost more than a linear limiter; shocks and admissibility scaling reduce the smooth formal order.',
'cweno5-char':'A wider stencil and eigenbasis work increase cost; limiting near jumps can remove the formal higher-order advantage.',
'weno5-z':'Adds a global smoothness correction but still requires a wider stencil and admissibility control; improved weights do not guarantee better discontinuous profiles.',
'teno3':'Binary selection depends on the detector and threshold. This educational variant has no blanket accuracy or positivity guarantee.',
'teno5':'Binary stencil rejection depends on cutoff and smoothness data; positivity control and adequate boundary treatment are still required.',
'teno7':'Wider stencils and binary cutoff decisions add sensitivity and work. This educational variant is not TENO-A/LAD.',
'muscl-thinc-bvd':'Single-stage componentwise selection can make imperfect interface choices and needs positivity controls. It is not the full published multistage BVD family.'}
def card(key):
    title,algorithm,limit=methods[key];seen.append(key)
    limit=LIMIT.get(key,limit)
    t=Table([[p(title,'card')],[p('<b>Algorithm.</b> '+algorithm)],[p('<b>Advantage.</b> '+ADV[key])],[p('<b>Limitation.</b> '+limit)]],colWidths=[WIDTH])
    t.setStyle(TableStyle([('BACKGROUND',(0,0),(-1,0),LIGHT),('LINEBEFORE',(0,0),(0,-1),2,TEAL),('LEFTPADDING',(0,0),(-1,-1),10),('RIGHTPADDING',(0,0),(-1,-1),10),('TOPPADDING',(0,0),(-1,0),7),('BOTTOMPADDING',(0,-1),(-1,-1),4)]))
    story.extend([KeepTogether([t]),Spacer(1,10)])
class CellsDiagram(Flowable):
    def __init__(self):Flowable.__init__(self);self.width=WIDTH;self.height=123
    def draw(self):
        c=self.canv;x0=35;cw=138;y=46
        for i,label in enumerate(['cell i-1','cell i','cell i+1']):
            c.setFillColor(LIGHT);c.setStrokeColor(TEAL);c.rect(x0+i*cw,y,cw,44,fill=1)
            c.setFillColor(BLUE);c.setFont('Lecture',10);c.drawCentredString(x0+(i+.5)*cw,y+25,label);c.drawCentredString(x0+(i+.5)*cw,y+9,'conserved mean')
        c.setFont('LectureBold',9)
        for x in [x0+cw,x0+2*cw]:
            c.drawCentredString(x,y-15,'shared face flux');c.setStrokeColor(TEAL);c.line(x-20,y+56,x+20,y+56);c.line(x+20,y+56,x+14,y+60);c.line(x+20,y+56,x+14,y+52)
class NumberedCanvas(canvas.Canvas):
    def __init__(self,*args,**kwargs):super().__init__(*args,**kwargs);self.saved=[]
    def showPage(self):self.saved.append(dict(self.__dict__));self._startPage()
    def save(self):
        n=len(self.saved)
        for state in self.saved:
            self.__dict__.update(state);self.setFillColor(GREY);self.setFont('Lecture',8)
            self.drawString(44,26,'Ehsan Roohi | CFD II | English teaching edition | 08 October 2026');self.drawRightString(W-44,26,f'{self._pageNumber} / {n}')
            super().showPage()
        super().save()
class WaveFans(Flowable):
    def __init__(self):Flowable.__init__(self);self.width=WIDTH;self.height=195
    def draw(self):
        c=self.canv;panel=WIDTH/3
        for j,title in enumerate(['Physical fan','HLL: two bounds','HLLC: add contact']):
            x=j*panel;origin=x+panel*.48;y=35;top=150
            c.setFillColor(BLUE);c.setFont('LectureBold',10);c.drawCentredString(x+panel/2,174,title)
            c.setStrokeColor(GREY);c.setLineWidth(.5);c.line(x+8,y,x+panel-6,y);c.line(origin,y,origin,top+8)
            c.setFont('Lecture',8);c.drawString(origin+4,top+5,'t');c.drawString(x+panel-10,y-10,'x')
            c.setStrokeColor(TEAL);c.setLineWidth(1.7)
            c.line(origin,y,x+20,top);c.line(origin,y,x+panel-12,top)
            if j==0:
                c.setLineWidth(.7)
                for z in [25,31,37,43,49]:c.line(origin,y,x+z,top)
            if j in (0,2):
                c.setStrokeColor(colors.HexColor('#B75E21'));c.setDash(4,2);c.line(origin,y,x+panel*.67,top);c.setDash()
            c.setFillColor(GREY);c.setFont('Lecture',8);c.drawCentredString(x+panel/2,10,'schematic wave speeds')
def header(c,doc):
    c.saveState();c.setStrokeColor(TEAL);c.setLineWidth(.6);c.line(44,H-37,W-44,H-37);c.setFillColor(GREY);c.setFont('Lecture',8);c.drawString(44,H-29,'MODERN CFD ALGORITHMS | FROM FACE FLUXES TO VERIFIABLE RESULTS');c.restoreState()

heading('CFD II / lecture notes','Modern CFD algorithms<br/>for the cone and shock tube')
add('A teaching companion to the 2026 Python notebooks','sub');story.append(Spacer(1,18))
box('The question this lecture answers','How do numerical methods move mass, momentum and energy across a cell face, and why do their computed shocks and contacts differ?')
add('Ehsan Roohi<br/>CFD II project modernization<br/>English edition - 08 October 2026');story.append(Spacer(1,14))
bullets(['15 implemented Euler face fluxes + standalone JST.','20 finite-volume reconstructions and DG degrees 0, 1 and 2.','Algorithm steps, advantages, limitations and coding exercises.','Actual Sod results, fixed uniform-grid refinement and exact-reference diagnostics.'])
box('Scope','This lecture teaches the algorithms selectable in the shared cone/shock-tube code. The plotted benchmark solves inviscid planar Euler equations. It does not validate physical viscous transport or cover every CFD method in use.')
add('Source baseline: public repository commit e7cd1655a07cbae51ab98d22bec87c274393fa55. Numerical algorithms and recorded results are unchanged by this lecture.','small')

heading('01 / teaching route','How to teach and use these notes')
table([['Session','Topics','Notebook activity'],['1: equations and updates','Euler waves, cell means, conservation, SSPRK3','Read Sod profiles; implement a local LF flux.'],['2: face solvers','Roe, HLL/HLLC, splitting, AUSM, entropy/JST','Hold CWENO3 fixed and compare fluxes.'],['3: high-order spatial methods','MUSCL, CWENO, WENO, TENO, BVD, DG','Hold HLLC fixed and compare reconstructions.'],['4: evidence and decisions','Admissibility, grid errors, interface width, failure labels','Compare six uniform grids and explain a method choice.']],[72,211,WIDTH-283])
sub('Learning outcomes')
bullets(['Trace a conservative update from reconstructed states to one shared flux.','Distinguish flux, reconstruction, time integration and physical models.','Explain why smooth formal order does not predict shock-profile accuracy.','Compare errors, positivity, conservation and wave resolution together.'])
sub('Reading route')
add('Foundations → face fluxes → reconstruction → DG and safeguards → measured comparisons → laboratory and exercises. Each named method has an algorithm, an advantage and a limitation.')
box('Prerequisites','Conservation laws, ideal-gas pressure, basic eigenvectors and Python arrays. The finite-volume update is developed here; no packaged CFD solver is required.')

heading('02 / governing equations','One physical system, many algorithms')
eq('U = (ρ, ρu, E)<super>T</super><br/>F(U) = (ρu, ρu<super>2</super> + p, u(E+p))<super>T</super><br/>∂U/∂t + ∂F/∂x = 0')
eq('p = (γ-1)[E - ρu<super>2</super>/2], &nbsp; a = √(γp/ρ), &nbsp; H = (E+p)/ρ')
bullets(['ρ is density; u is normal velocity; E is total energy per volume.','Pressure closes the ideal-gas equations; a sets an acoustic-speed scale.','Euler wave speeds are u-a, u and u+a in planar flow.','The shared code carries transverse momentum too; it is zero in Sod.'])
box('Physical versus numerical dissipation','The Euler comparisons have no molecular-viscosity term. Damping from LF, JST or a limiter is numerical dissipation. A new flux does not add physical viscosity or a turbulence model.')
sub('Sod initial data')
table([['Region','ρ','u','p'],['Left of x=0.5','1','0','1'],['Right of x=0.5','0.125','0','0.1']],[200,90,90,WIDTH-380])
add('Domain [0,1]; γ=1.4; requested time t=0.2; outflow boundaries.','small')

heading('03 / finite volume','Build a conservative cell update')
story.append(CellsDiagram());eq('dŪ<sub>i</sub>/dt = -(Fhat<sub>i+1/2</sub> - Fhat<sub>i-1/2</sub>)/Δx')
bullets(['Integrate the conservation law over a cell and store its conserved mean.','Reconstruction supplies two states at each face.','The face solver turns these two states into one flux.','Use the identical flux in both neighbors. Interior contributions cancel in the domain budget.'])
sub('A small mass-balance calculation')
add('Let Δx=0.01, Δt=0.001, left mass flux=0.20 and right mass flux=0.25. The Euler substep changes density by -(0.001/0.01)(0.25-0.20)=-0.005. A density of 1 becomes 0.995. After summing cells, only boundary flux contributions remain.')
box('Coding check','Compute each face flux once into an array. Check a two-cell update by hand. Updating primitive variables independently generally does not reproduce the conservative Euler update.')

heading('04 / wave physics','What must a shock tube resolve?')
table([['Feature','Physical signature','Numerical challenge'],['Rarefaction','Continuous expansion; density and pressure vary smoothly.','Retain smooth accuracy without flattening the fan.'],['Contact','Density jump; pressure and normal velocity continuous.','Limit spreading of a transported density interface.'],['Shock','Compression jump in density, pressure and velocity.','Capture strength and position without ringing.']],[100,195,WIDTH-295])
sub('Why use cone flow and shock tubes?')
add('Cone flow couples wave capturing to geometry, source integration, boundaries and a structured spatial flow. A planar shock tube removes the cone source and isolates transient wave interactions. These are complementary tests that stress different components.')
box('Interpret a rounded transition','A finite-volume jump spans finite cells. Refinement should reduce its physical thickness. A rarefaction is physically smooth and should not be sharpened into a step.')
bullets(['Sod: expansion, contact and shock.','Lax: stronger waves.','Double rarefaction: low-pressure/density stress.','Stationary contact: contact diffusion without transport.'])

heading('04b / visual reasoning','From a wave fan to a face flux')
story.append(WaveFans())
add('At a local face, left/right states define a Riemann problem. Rays have slope x/t equal to a wave speed. The flux samples the state at the face ray x/t=0. This original schematic shows topology, not the measured Sod wave speeds.')
table([['Description','Retained information'],['Exact local solver','Shock/rarefaction branches and a contact; nonlinear fan structure.'],['HLL approximation','Two outer bounds and one intermediate state; contact not represented separately.'],['HLLC approximation','Outer bounds plus a contact dividing two star states.']],[165,WIDTH-165])
box('Classroom prompt','Where is x/t=0 in each panel? If the entire fan travels right, which physical state supplies the face flux? Answer: the left state. Repeat with an entirely left-going fan.')
add('Teaching idea: use local wave diagrams before full numerical profiles, following Ketcheson, LeVeque and del Razo, Riemann Problems and Jupyter Solutions, Euler approximate-solvers chapter. See the linked reading map at the end.','small')

heading('05 / time integration','Advance with SSPRK3')
eq('A = U<super>n</super> + Δt L(U<super>n</super>)<br/>B = 3U<super>n</super>/4 + [A + Δt L(A)]/4<br/>U<super>n+1</super> = U<super>n</super>/3 + 2[B + Δt L(B)]/3')
eq('Δt<sub>FV</sub> = CFL · Δx / max(|u|+a)')
bullets(['Recompute the full spatial operator at each stage, including boundaries.','Clip Δt to reach the requested time. Recorded uniform-grid CFL is 0.25.','Reject invalid stage means and retry at a smaller step.','Strong-stability-preserving terminology does not guarantee positivity or no oscillations for every spatial method.'])
box('Worked estimate','For Δx=0.01, max(|u|+a)=2 and CFL=0.25, Δt=0.00125. Halving Δx halves the allowable step. Refinement increases both cell work and step count.')
add('Step pseudocode: choose Δt → reconstruct admissible states → evaluate fluxes → SSPRK stage → validate means → repeat → accept, or restore the old state and halve Δt. Retry counts are finite.')

heading('06 / Riemann solvers','Exact and linearized local wave fans')
for k in ['godunov','roe','roe-nc']:card(k)
box('Roe dissipation','Fhat = (F<sub>L</sub>+F<sub>R</sub>)/2 - ½ Σ |λ<sub>k</sub>| α<sub>k</sub> r<sub>k</sub><br/>Decompose the jump into waves r<sub>k</sub>; their speed magnitudes set wave-dependent dissipation.')

heading('07 / HLL family','Two waves, then restore the contact')
for k in ['hlle','hllc']:card(k)
eq('S<sub>L</sub> = min(u<sub>L</sub>-a<sub>L</sub>, u<sub>R</sub>-a<sub>R</sub>, 0)<br/>S<sub>R</sub> = max(u<sub>L</sub>+a<sub>L</sub>, u<sub>R</sub>+a<sub>R</sub>, 0)<br/>Fhat<sub>HLL</sub> = [S<sub>R</sub>F<sub>L</sub>-S<sub>L</sub>F<sub>R</sub>+S<sub>L</sub>S<sub>R</sub>(U<sub>R</sub>-U<sub>L</sub>)]/(S<sub>R</sub>-S<sub>L</sub>)')
box('Name-to-code caution','The hlle key uses Davis min/max acoustic bounds, not Einfeldt/Roe-averaged bounds. The label alone does not specify the implemented wave-speed estimate.')

heading('08 / simpler upwinding','LF and flux-vector splitting')
for k in ['rusanov','global-lf','vanleer','steger-warming']:card(k)

heading('09 / AUSM foundation','Separate mass transport and pressure')
eq('Fhat = ṁ · (1, u<sub>up</sub>, H<sub>up</sub>)<super>T</super> + (0, p<sub>face</sub>, 0)<super>T</super>')
add('The mass-flux sign chooses upwind velocity and enthalpy. Split interface pressure acts on normal momentum. The four-component code also transports transverse momentum.')
for k in ['ausm','ausm-plus']:card(k)
box('What changes within the family?','Interface sound speed, Mach/pressure polynomials and dissipation corrections. Conserved variables, equation of state and physical Euler flux remain the same.')

heading('10 / all-speed splitting','AUSM+-up, AUSM+-up2 and SLAU2')
for k in ['ausm-up','ausm-up2','slau2']:card(k)
box('Mechanism before ranking','Pressure-jump corrections change mass/pressure coupling; velocity terms change pressure dissipation. All-speed intent does not imply the smallest Sod error. Hold reference-Mach settings fixed in grid comparisons.')

heading('11 / central constructions','Entropy core and JST dissipation')
card('ec-lf')
box('jst - standalone central stencil','<b>Algorithm.</b> Evaluate F at the averaged face state; subtract speed-scaled second/fourth-difference dissipation. Pressure curvature activates stronger second-difference damping near shocks.<br/><br/><b>Advantage.</b> The sensor exposes the damping mechanism explicitly.<br/><br/><b>Limitation.</b> The current stencil can ring near jumps. It is numerical damping, not physical viscosity, and is not interchangeable with DG or every reconstruction.')
eq('ε<sub>2</sub> = 0.5 s<sub>p</sub>, &nbsp; ε<sub>4</sub> = max(0, 0.02-ε<sub>2</sub>)')
box('Different entropy claims','A central face identity, dissipative two-state flux, reconstructed semidiscrete scheme and fully discrete boundary-treated scheme are distinct claims. Face checks do not prove global entropy stability.')

heading('12 / high-order reconstruction','Cell means are not nodal samples')
eq('Ū<sub>i</sub> = (1/Δx) ∫<sub>cell i</sub> P<sub>i</sub>(x) dx<br/>U<sub>i+1/2</sub><super>L</super> = P<sub>i</sub>(x<sub>i+1/2</sub>)')
bullets(['Fit moments to cell means; evaluate the fitted representation at faces.','A limiter or nonlinear weight suppresses candidates that cross a jump poorly.','Higher formal order requires smooth information; a wider stencil does not make a discontinuity smooth.','Characteristic reconstruction projects into local wave variables, reconstructs and transforms back.'])
box('Three independent controls','Flux: interface dissipation. Reconstruction: face data. Mesh spacing: physical resolution. Changing all three at once obscures why an answer changes.')
sub('A fair experiment')
add('Hold HLLC, grid, boundaries, final time and SSPRK3 fixed. Change only reconstruction. Inspect contact width, overshoot and density/pressure errors, then repeat matched grids before claiming observed order.')

heading('13 / linear baselines','First order, MC and minmod')
for k in ['first','muscl','muscl-minmod']:card(k)
sub('Hand calculation')
add('Means 1.0, 1.2, 1.3 give a=0.2, b=0.1. MC slope=min(2a,2b,(a+b)/2)=0.15; minmod=0.1. On normalized coordinates [-½,½], MC face states are 1.125 and 1.275; minmod gives 1.15 and 1.25. Both preserve mean 1.2.')

heading('14 / MUSCL choices','Smoothness versus compression')
for k in ['muscl-vanleer','muscl-superbee','muscl-vanalbada']:card(k)
box('Same data, different slopes','For a=0.2, b=0.1: Van Leer=0.13333; Superbee=0.2; Van Albada≈0.12. When signs disagree, these implemented variants return zero. Sharper permitted slopes can be more compressive.')

heading('15 / candidates','ENO and cell-based CWENO')
for k in ['eno2','cweno3','cweno5']:card(k)
eq('ω<sub>k</sub> = α<sub>k</sub>/Σα<sub>j</sub>, &nbsp; α<sub>k</sub> = d<sub>k</sub>/(β<sub>k</sub>+ε)<super>2</super>')
add('d are optimal weights; β measure derivative smoothness. A central CWENO candidate allows the linear combination to recover the optimal wide-stencil polynomial.','small')

heading('16 / characteristics','Weight waves instead of components')
for k in ['cweno3-char','cweno5-char']:card(k)
bullets(['Evaluate a local Euler eigenbasis at the cell state.','Project candidate coefficients into wave variables.','Compute nonlinear weights for each wave component.','Transform back and apply admissibility controls.'])
box('Advantage and cost','Separate wave weighting can reduce acoustic/contact coupling. It requires eigenvectors and extra transformations. It alone guarantees neither positive pressure nor no overshoot.')

heading('17 / WENO weighting','Jiang-Shu versus Z weights')
for k in ['weno3-js','weno3-z','weno5-js','weno5-z']:card(k)

heading('17b / worked weights','Why a troubled stencil loses influence')
add('This synthetic arithmetic example illustrates JS weighting. It is not a fitted Euler stencil or an accuracy experiment. Let candidate face values be q=(1,3,10), optimal weights d=(0.1,0.6,0.3), indicators β=(0.01,0.02,1) and ε=10<super>-6</super>.')
eq('α<sub>k</sub> = d<sub>k</sub>/(β<sub>k</sub>+ε)<super>2</super>, &nbsp; ω<sub>k</sub> = α<sub>k</sub>/Σα<sub>j</sub>')
ds=[.1,.6,.3];bs=[.01,.02,1];qs=[1,3,10]
alphas=[d/(b+1e-6)**2 for d,b in zip(ds,bs)];weights=[a/sum(alphas) for a in alphas]
table([['Candidate','q','β','d','α','ω']]+[[k,qs[k],bs[k],ds[k],f'{alphas[k]:.5f}',f'{weights[k]:.7f}'] for k in range(3)],[65,55,65,55,125,WIDTH-365])
add(f'Linear combination: Σd q = {sum(d*q for d,q in zip(ds,qs)):.4f}.<br/>Nonlinear combination: Σω q = {sum(w*q for w,q in zip(weights,qs)):.6f}. The large indicator on candidate 2 makes its weight very small.')
box('Smooth-region limit','When all candidate smoothness indicators are approximately equal, the normalized JS weights approach the optimal linear weights d. When a candidate crosses a jump, the indicator should suppress it. Near critical points, the route back to optimal weights can be imperfect; Z weighting addresses smooth-weight recovery differently.')
add('Teaching idea: explain reconstruction, smoothness indicators and normalization as separate operations, following Shu, NASA/CR-97-206253, Sections 2.1-2.3. The numerical example is newly constructed for this lecture; actual code ε settings remain as documented.','small')

heading('18 / wider stencils','WENO7 and targeted TENO')
for k in ['weno7-js','teno3','teno5','teno7']:card(k)
add('TENO cutoff=10<super>-5</super>; exponent=6. TENO3/7 are educational inverse-smoothness cutoff variants. No TENO-A/LAD or WENO7-Z option is implemented.','small')

heading('19 / interface candidates','THINC and boundary variation diminishing')
card('muscl-thinc-bvd')
bullets(['Build MC-linear and bounded tanh (THINC) face pairs.','Apply the implemented monotone-component conditions.','Compare boundary-jump variation with neighbor candidates.','Choose the lower variation componentwise, then control face admissibility.'])
box('Decision exercise','Linear candidate boundary jumps 0.08 and 0.05 sum to 0.13. THINC jumps 0.03 and 0.04 sum to 0.07. This local comparison selects THINC. The actual code considers neighbor candidates; the arithmetic illustrates the principle.')
box('Roles remain separate','BVD changes face states; the selected Riemann flux still computes transport. This is a single-stage educational selection, not a complete published multistage P4-THINC-BVD implementation.')

heading('20 / DG','Evolve a polynomial in each element')
eq('U<sub>h</sub>(ξ,t) = Σ<sub>k=0…p</sub> Û<sub>k</sub>(t) P<sub>k</sub>(ξ), &nbsp; ξ in [-1,1]')
add('Multiply by basis functions and integrate by parts. Volume quadrature handles physical flux; common numerical fluxes couple element boundaries. Legendre coefficients make modal degrees of freedom explicit.')
eq('∫<sub>K</sub> (∂U<sub>h</sub>/∂t) v dx = ∫<sub>K</sub> F(U<sub>h</sub>) v′ dx - [Fhat v]<sub>∂K</sub>')
table([['Degree','Stored content','Smooth order / limitation'],['0','Cell mean','First order; matches first-order FV for matched steps and flux.'],['1','Mean + linear mode','Second order; limiting reduces order near shocks.'],['2','Mean + linear + quadratic modes','Third order; extra work and smaller stable time steps.']],[55,170,WIDTH-225])
bullets(['All 15 pointwise fluxes available; JST excluded.','At least max(3,p+2) Gauss points; CFL restriction includes 1/(2p+1).','Characteristic TVB controls slopes; altered slopes zero higher modes. Recorded TVB parameter=0.','Sampled positivity scaling reduces higher modes toward the conserved mean.'])
box('Consistency check','Degree-0 DG agrees with first-order FV to roundoff for all shared fluxes. This checks their matched limit; it does not establish higher-degree shock accuracy.')
add('Teaching connection: set test function v=1, so v′=0. The volume term vanishes and the mean evolves from boundary fluxes, recovering FV. This route follows Persson’s Math 228B DG lectures; our code uses a modal Legendre representation.','small')

heading('21 / safeguards','Positivity and conservation are different')
bullets(['Positive density is insufficient: internal energy after kinetic subtraction must give positive pressure.','Scale toward an admissible mean instead of independently clipping conserved components.','Reducing DG higher modes preserves the mean.','An invalid stage mean requires rejecting/reducing the step; reconstruction scaling cannot fix it.'])
eq('U<sub>limited</sub>(x) = Ū + θ[U<sub>candidate</sub>(x)-Ū], &nbsp; 0 ≤ θ ≤ 1')
box('Sampled controls','Polynomial/face states are controlled at sample points. Stage means are tested with bounded retry counts. This does not prove admissibility at every point of every polynomial.')
sub('Label failures honestly')
add('The original study has 592 runs: 589 completed integrity-passing runs and three earlier-time double-rarefaction failures. Preserve actual stopping times. Do not assign final-time errors to early states.')

heading('22 / flow properties','Read actual computed fields')
fig('uniform-physical-properties.png','Executed HLLC + CWENO3 at t=0.2 on six uniform grids. Curves connect computed cell means; the dashed line is the exact reference.',430)
bullets(['Density reveals contact and shock; pressure distinguishes the shock from the contact.','Derived internal energy, Mach number and temperature come from the same conserved states.','Keep the physical rarefaction smooth as discontinuous transitions narrow.'])
add('With nondimensional R=1: T=p/ρ; specific internal energy=p/[(γ-1)ρ]; Mach=|u|/a.','small')

heading('23 / mesh refinement','Change physical resolution uniformly')
fig('uniform-shock-contact-detail.png','HLLC + CWENO3. N=80,160,320,640,1280,2560; CFL=0.25, t=0.2 fixed. No adaptation or cosmetic sharpening.',355)
table([['N','Density L1','Shock width','Contact width'],['80','0.00881958','0.0355518','0.0600959'],['320','0.00239207','0.00906783','0.0225290'],['1280','0.000694095','0.00225816','0.00807732'],['2560','0.000416369','0.00112936','0.00481002']],[55,145,145,WIDTH-345])
add('HLLC density L1 falls 95.28%. Shock width shrinks about 31.5-fold; contact width 12.5-fold. The finest shock spans 2.89 cells and contact 12.31 cells. Physical narrowing need not mean fewer transition cells.')

heading('24 / all methods','Compare waves, not only one error')
fig('uniform-all-method-shock-detail.png','All 15 pointwise fluxes use CWENO3; JST uses its own stencil. Six uniform grids in every panel, with the exact Sod shock.',485)
box('An instructive exception','The current JST stencil shows oscillations at the jump. A narrow 10-90% width can miss this weakness. Compare width, overshoot, field errors and admissibility together; preserve unfavorable results too.')

heading('25 / convergence','Formal order and measured error differ')
fig('uniform-grid-error-convergence.png','Actual density, velocity and pressure L1 errors; each algorithm is held fixed as Δx decreases.',330)
eq('L1(q) = Σ<sub>i</sub> |q<sub>i,num</sub>-q<sub>i,ref</sub>| Δx<br/>Observed order r = log(e<sub>N</sub>/e<sub>2N</sub>)/log(2)')
bullets(['Integrate exact conserved states into cell means, then convert to primitive variables for the reference.','Discontinuities lower global observed order even for high-order smooth reconstruction.','Separate mesh-error convergence from transient RHS/time histories; this is not a steady residual solve.'])
box('Actual coverage','96 comparisons=16 methods × 6 grids, all integrity-passing. The original full FV matrix=15 fluxes × 20 reconstructions=300 runs on 80 cells. Six-grid refinement does not repeat all 300 combinations or DG.')

heading('26 / laboratory','Choose methods through controlled tests')
table([['Question','Hold fixed','Change / inspect'],['Flux effect on contact','CWENO3, grid, time, CFL','HLL/HLLC/LF; density width and errors.'],['Reconstruction effect','HLLC, grid, time, CFL','First/MUSCL/WENO; errors and overshoot.'],['Steepening benefit','HLLC, stepper, grid','MC/BVD; width plus errors and failures.'],['Higher DG degree','Flux and final time','Degree, limiting, time step and error.'],['Enough resolution?','Algorithm and parameters','Uniform N; physical widths and errors.']],[150,135,WIDTH-285])
sub('Student deliverable')
bullets(['Six-property plot, exact reference and method/grid legend.','Shock/contact zoom and L1/transition-width table.','Uniform-grid errors and observed orders with clearly stated norms.','Admissibility, conservation and failure table; explain one favorable and one unfavorable result.'])
box('Teach a simple baseline first','First-order FV + Rusanov is short enough to debug. HLLC + CWENO3 is a useful controlled comparison here. Add one change at a time; no method is declared universally best.')

heading('27 / questions','Exercises before reading the answers')
bullets(['1. Sum mass updates over three cells. Which flux terms remain?','2. CFL=0.25, Δx=1/160, max(|u|+a)=2: find Δt. What changes at N=320?','3. With a=0.2, b=0.1, find MC, minmod, Van Leer, Superbee and Van Albada slopes. Find faces for mean=1.2.','4. Why can HLL diffuse a contact more than HLLC with identical reconstruction?','5. Why does an exact Godunov face flux not yield an exact full shock-tube computation?','6. Density L1 is 0.00881958 at N=80 and 0.00453499 at N=160. Estimate observed order.','7. Can an oscillatory shock have a narrower 10-90% width and still be worse? Name two additional checks.','8. Design a fair WENO-JS versus WENO-Z experiment.','9. Why do seventh-order faces not establish seventh-order cone accuracy with a CWENO5 source polynomial?','10. Why does degree-0 DG agreement with FV not validate degree-2 shock accuracy?'])

heading('28 / answers','Answers and discussion prompts')
bullets(['1. Shared interior terms cancel; only boundary fluxes remain. SSPRK conservation uses stage-weighted boundary contributions.','2. Δt=0.25/(160×2)=0.00078125. At N=320: 0.000390625.','3. MC=0.15; minmod=0.10; Van Leer≈0.13333; Superbee=0.20; Van Albada≈0.12. Face values=1.2±s/2.','4. HLL retains outer bounds only. HLLC adds the contact speed and star states.','5. Reconstruction, finite cells, time integration, limiting and boundaries still introduce error. Exact refers only to the local Riemann problem.','6. log(0.00881958/0.00453499)/log(2)≈0.96. This is not CWENO3 smooth formal order.','7. Yes. Check overshoot/undershoot and exact-reference errors, as well as positivity/conservation.','8. Hold flux, grid, CFL, time, boundaries and characteristic treatment fixed; change only weights. Then repeat matched grids.','9. Source integration participates in the whole spatial operator. Higher face order cannot override a lower-order source component.','10. Degree 0 checks the mean update. Degrees 1/2 additionally require volume quadrature, modal evolution, limiting and appropriate time steps.'])

heading('29 / sources','Connect lessons to runnable code')
url='https://github.com/Ehsan-Roohi/Ehsan-Roohi/tree/main/courses/cfd-2/modern-python'
add(f'<link href="{url}" color="#087F8C">Open the course code and teaching files</link>')
table([['Topic','Shared conical/ module'],['FV / DG / cases','shock_suite.py and shared DG core'],['Roe, HLL/C, LF, Van Leer','flux.py'],['Godunov, uncorrected Roe, SW, JST','additional_fluxes.py'],['AUSM / SLAU2','ausm.py'],['Entropy core + LF','entropy_flux.py'],['CWENO / MUSCL / characteristic weights','reconstruction.py'],['ENO / WENO / TENO / BVD','advanced_reconstruction.py'],['Recorded grids and diagnostics','results-uniform/catalog.json; shock-tube/UNIFORM_GRID_RESULTS.md']],[215,WIDTH-215])
sub('Primary reading and code comparisons')
for title,link in [('NASA shocktube: AUSM/AUSM+ examples','https://github.com/nasa/shocktube'),('Sun, Inaba and Xiao: Boundary Variation Diminishing (2016)','https://arxiv.org/abs/1602.00814'),('Ihme Group Quail: DG teaching and prototyping','https://github.com/IhmeGroup/quail'),('Riemann-Solvers: exact and approximate examples','https://github.com/cangyu/Riemann-Solvers')]:
    add(f'<link href="{link}" color="#087F8C">{html.escape(title)}</link>','small')
add('Descriptions and numerical values follow this course implementation and measured records. External codes/papers are reading references; their full algorithms and guarantees are not automatically implemented here.','small')
box('Report with every numerical claim','Equations, flux, reconstruction, grid, boundaries, stepper, CFL, parameters, actual stopping time, reference definition and integrity checks. Preserve attribution and distinguish formal smooth order from measured error.')

heading('30 / open lecture reading map','What inspired the teaching structure?')
readings=[
('LeVeque / University of Washington','Finite Volume Methods for Hyperbolic Problems: 2023 supplementary lectures and PDF slides','https://www.clawpack.org/fvmhp_materials/','Teaching route: conservation and waves before high-resolution comparisons; pair explanations with runnable experiments.'),
('Ketcheson, LeVeque and del Razo','Riemann Problems and Jupyter Solutions: Euler approximate solvers','https://www.clawpack.org/riemann_book/html/Euler_approximate.html','Teaching route: local wave fans first, then approximate solvers and their failure mechanisms. The schematic here is newly drawn.'),
('Shu / Brown; NASA ICASE report (1997)','ENO and WENO lecture notes, Sections 2.1-2.3','https://academicweb.nd.edu/~zxu2/acms60790S13/Shu-WENO-notes.pdf','Teaching route: cell-average reconstruction, candidate stencils and nonlinear weights. The worked weight arithmetic here is new.'),
('Persson / UC Berkeley Math 228B','Discontinuous Galerkin Methods for Conservation Laws','https://persson.berkeley.edu/math228b/slides/dg_slides.pdf','Teaching route: derive FV from piecewise constants, then extend to element polynomials and common boundary fluxes.')]
for author,title,link,idea in readings:
    sub(author);add(f'<link href="{link}" color="#087F8C">{html.escape(title)}</link>','small');add(idea)
add('Access checked 08 October 2026. These notes use original wording, exercises, diagrams and the course’s own measured figures. External lecture PDFs are linked, not reproduced. Their finite-difference, adaptive or viscous extensions are not added to this benchmark.','small')
add('<b>Michigan context:</b> Fidkowski’s official teaching page lists CFD I/II, but that page does not supply a downloadable 2025/2026 lecture set. No such set is claimed as the source of these notes. <link href="https://public.websites.umich.edu/~kfid/teaching.html" color="#087F8C">Official teaching page</link>.','small')

assert set(seen)==set(methods) and len(seen)==35
assert len(plan)==33
doc=SimpleDocTemplate(str(PDF),pagesize=A4,leftMargin=44,rightMargin=44,topMargin=58,bottomMargin=48,title='Modern CFD Algorithms - English Lecture Notes',author='Ehsan Roohi',subject='Finite-volume fluxes, reconstruction, DG and measured shock-tube comparisons')
doc.build(story,onFirstPage=header,onLaterPages=header,canvasmaker=NumberedCanvas)
reader=PdfReader(PDF);texts=[page.extract_text() or '' for page in reader.pages]
assert len(texts)==len(plan),(len(texts),len(plan),'A teaching section overflowed')
full='\n'.join(texts);assert all(key in full for key in methods)
record={'pages':len(texts),'language':'English','fluxes':15,'standalone_jst':True,'fv_reconstructions':20,'dg_degrees':[0,1,2],'all_selectable_methods_present':True,'sha256':hashlib.sha256(PDF.read_bytes()).hexdigest(),'bytes':PDF.stat().st_size,'planned_sections':plan,'numerical_source_commit':'e7cd1655a07cbae51ab98d22bec87c274393fa55','created_date':'2026-10-08','simulations_rerun':False}
(OUT/'lecture-validation.json').write_text(json.dumps(record,indent=2)+'\n')
print(json.dumps(record,indent=2))
