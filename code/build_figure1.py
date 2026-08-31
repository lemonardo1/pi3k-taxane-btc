def build_fig1():
    fig = newfig(7.2, 8.0)
    gsA = fig.add_gridspec(3, 2, height_ratios=[1.45, 2.25, 2.30],
                           left=0.105, right=0.975, top=0.945, bottom=0.065,
                           hspace=0.62, wspace=0.34)
    # --- (a) schematic
    axa = fig.add_subplot(gsA[0, :]); axa.set_axis_off(); axa.set_xlim(0,100); axa.set_ylim(0,100)
    def box(x,y,w,h,txt,fc,ec=None,fs=6.4):
        axa.add_patch(FancyBboxPatch((x,y),w,h,boxstyle="round,pad=0.6,rounding_size=2",
                     lw=0.7, fc=fc, ec=ec or fc, zorder=1))
        axa.text(x+w/2, y+h/2, txt, ha='center', va='center', fontsize=fs, color=INK, zorder=2, linespacing=1.35)
    def arr(x,y0,y1):
        axa.annotate('', xy=(x,y1), xytext=(x,y0),
                     arrowprops=dict(arrowstyle='-|>', lw=0.8, color=INK, mutation_scale=7))
    box(1.5, 70, 29.5, 27, "Index cohort (Yonsei)\n$n$ = 287 · GAP 198 / GC 89\nTruSight Oncology 500", '#F6DDE0', '#D8A9B0')
    box(34.5,70, 29.5, 27, "Validation cohort (CHA)\n$n$ = 167 · GAP 103 / GC 64\nOncomine panels", '#F6DDE0', '#D8A9B0')
    box(68,  70, 30.5, 27, "Public BTC cohorts ($n$ = 10)\n1,356 profiled · mutation calls\n560 with taxane exposure known", '#D9E4F2', '#A8BFDC')
    box(1.5, 36, 62.5, 20, "Treatment-effect modification\npathway screen · IPTW · T-learner", '#EFEFEF', '#CFCFCF')
    box(68,  36, 30.5, 20, "Survival association\n(taxane exposure documented)", '#EFEFEF', '#CFCFCF')
    for x in (16.2, 49.2, 83.2): arr(x, 69, 57.5)
    box(12, 2, 76, 21, "ESM-2 variant scoring → mode-aware continuous PI3K activation score\n"
                       "56 observed variants · 283,442-variant in-silico background (14 genes)",
        '#E7E0F0', '#BFB2D4')
    arr(33, 35, 24); arr(83.2, 35, 24)
    axa.set_title('Study design: three data sources, one continuous score', loc='left', pad=4)
    L(fig, axa, 'a', dx=-0.085, dy=-0.004)

    # --- (b) prevalence
    axb = fig.add_subplot(gsA[1, 0])
    prev = [("Public, mutation\ncalls\n$n$=1,356", 220/1356*100, BLUE),
            ("Public, taxane-\nunexposed\n$n$=560",  92/560*100, BLUE),
            ("Index\n(Yonsei)\n$n$=287",            58/287*100, RED),
            ("Validation\n(CHA)\n$n$=167",          52/167*100, RED)]
    xs = np.arange(4)
    axb.bar(xs, [p[1] for p in prev], color=[p[2] for p in prev], width=0.68, zorder=3)
    for x,(lab,val,c) in zip(xs, prev):
        axb.text(x, val+0.9, f"{val:.1f}%", ha='center', va='bottom', fontsize=7, color=INK)
    axb.set_xticks(xs); axb.set_xticklabels([p[0] for p in prev], fontsize=5.8)
    axb.set_ylabel('PI3K pathway alterations (%)')
    axb.set_ylim(0, 38); axb.margins(x=0.06)
    axb.set_title('Mutation calls give a public prevalence\nconcordant with the study cohorts', loc='left', pad=5)
    axb.axhline(220/1356*100, ls=':', lw=0.7, color=GREY, zorder=2)
    axb.text(3.42, 220/1356*100+0.8, 'public mean 16.2%', fontsize=6, color=GREY, va='bottom', ha='right')
    L(fig, axb, 'b')

    # --- (c) pathway screen
    ps = D['pathway_interaction_screen']
    pos = ps[ps.endpoint=='OS'].copy().sort_values('HR_int')
    namemap={'PI3K':'PI3K','HIPPO':'Hippo','WNT':'WNT','Cell_cycle':'Cell cycle','TGFB':'TGF-β',
             'TP53':'TP53','NOTCH':'Notch','RTK_RAS':'RTK–RAS','NRF2':'NRF2'}
    axc = fig.add_subplot(gsA[1, 1])
    n=len(pos)
    axc.vlines(1, -0.55, n-0.45, ls='--', lw=0.8, color=INK, zorder=1)
    for i,(_,r) in enumerate(pos.iterrows()):
        foc = r.pathway=='PI3K'; c = RED if foc else GREY
        axc.plot([r.lo, r.hi],[i,i], color=c, lw=1.9 if foc else 1.0, solid_capstyle='round', zorder=3 if foc else 2)
        axc.plot([r.HR_int],[i], 'o', ms=5.5 if foc else 3.8, color=c, zorder=4, mec='white', mew=0.6)
    axc.set_yticks(np.arange(n)); axc.set_yticklabels([namemap[p] for p in pos.pathway], fontsize=6)
    for tl, p in zip(axc.get_yticklabels(), pos.pathway):
        if p=='PI3K': tl.set_color(RED); tl.set_fontweight('bold')
    axc.set_xscale('log'); axc.set_xticks([0.25,1,4]); axc.set_xticklabels(['0.25','1','4'])
    axc.xaxis.set_minor_locator(mpl.ticker.NullLocator())
    axc.set_xlim(0.09, 11.0); axc.set_ylim(-2.05, n-0.45)
    axc.spines['left'].set_bounds(-0.45, n-0.45)
    axc.set_xlabel('Treatment × pathway interaction HR (OS)')
    axc.set_title('PI3K is the only pathway that\nmodifies GAP benefit', loc='left', pad=5)
    axc.text(0.098, -1.10, 'PI3K  HR 2.71 (1.40–5.27), $P$ = 0.0032', fontsize=6, color=RED, va='center', ha='left')
    axc.text(0.098, -1.72, '9 testable pathways, $n$ = 287', fontsize=6, color=GREY, va='center', ha='left')
    L(fig, axc, 'c')

    # --- (d) forest, own gridspec with wide left margin
    gsD = fig.add_gridspec(1, 1, left=0.235, right=0.975, top=0.318, bottom=0.075)
    fp = pd.read_csv('forest_plot_data.csv')
    axd = fig.add_subplot(gsD[0, 0])
    rows = [('hdr','Public cohorts, taxane exposure documented',None),
            ('dat','ICGC 2017 (Thai/European, surgical)',0),
            ('dat','MSK 2018 (first-line regimen names)',1),
            ('dat','TCGA (surgical)',2),
            ('dat','SMMU 2014 (Chinese, surgical)',3),
            ('dat','Pooled, cohort-stratified',4),
            ('dat','Random-effects meta-analysis',5),
            ('dat','+ chemotherapy-flag cohort (sensitivity)',6),
            ('hdr','By treatment arm',None),
            ('dat','MSK 2018 GC arm (no taxane)',7),
            ('dat','Yonsei GC arm (no taxane)',8),
            ('dat','Yonsei GAP arm (taxane-exposed)',9)]
    yv = np.arange(len(rows))[::-1]
    axd.vlines(1, -0.45, len(rows)-0.45, ls='--', lw=0.8, color=INK, zorder=1)
    for y,(kind,lab,idx) in zip(yv, rows):
        if kind=='hdr':
            axd.text(0.315, y, lab, fontsize=6.4, style='italic', color=INK, va='center', ha='left'); continue
        r = fp.iloc[idx]
        foc = 'GAP arm' in lab
        c = RED if foc else (BLUE if ('meta' in lab or 'GC arm' in lab) else GREY)
        mk = 'D' if ('meta' in lab or 'Pooled' in lab or 'sensitivity' in lab) else 'o'
        axd.plot([r.lo, r.hi],[y,y], color=c, lw=1.8 if foc else 1.1, solid_capstyle='round', zorder=3)
        axd.plot([r.hr],[y], mk, ms=6 if foc else 4.6, color=c, mec='white', mew=0.6, zorder=4)
        axd.text(6.4, y, f"{r.hr:.2f} ({r.lo:.2f}–{r.hi:.2f})", fontsize=6, va='center', ha='left',
                 color=c if foc else INK)
        axd.text(15.5, y, f"{int(r.n)} / {int(r.ev)}", fontsize=6, va='center', ha='right', color=GREY)
    axd.set_yticks(yv); axd.set_yticklabels([('' if k=='hdr' else l) for k,l,_ in rows], fontsize=6)
    axd.set_xscale('log'); axd.set_xlim(0.30, 16.0); axd.set_ylim(-1.35, len(rows)-0.45)
    axd.set_xticks([0.5,1,2,4]); axd.set_xticklabels(['0.5','1','2','4'])
    axd.xaxis.set_minor_locator(mpl.ticker.NullLocator())
    axd.set_xlabel('PI3K effect on overall survival (HR)')
    axd.text(6.4, len(rows)-0.75, 'HR (95% CI)', fontsize=6, style='italic', color=GREY, ha='left', va='center')
    axd.text(15.5, len(rows)-0.75, '$n$ / events', fontsize=6, style='italic', color=GREY, ha='right', va='center')
    axd.spines['bottom'].set_bounds(0.30, 4.6)
    axd.spines['left'].set_bounds(-0.45, len(rows)-0.45)
    axd.set_title('The PI3K association is stronger where a taxane was given than where it was not',
                  loc='left', pad=5)
    axd.text(0.315, -1.08, 'Ratio of taxane-exposed to taxane-unexposed HR = 1.74 (95% CI 1.08–2.80), '
             '$P$ = 0.023   ·   $I^2$ = 22%, $k$ = 4', fontsize=6, color=GREY, ha='left', va='center')
    L(fig, axd, 'd', dx=-0.215)
    return fig

fig1 = build_fig1()
fig1.savefig('fig/Figure1.png', dpi=350, facecolor='white')
ov,out = check(fig1,'Figure1')
