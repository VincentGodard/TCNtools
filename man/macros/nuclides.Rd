% Nuclide symbol with its mass number as a superscript: \nuclide{10}{Be}.
% <sup> in HTML (help pages and pkgdown site), a math superscript in the
% PDF manual, plain "10Be" in text help.
\newcommand{\nuclide}{\ifelse{html}{\out{<sup>#1</sup>#2}}{\ifelse{latex}{\eqn{^{#1}}#2}{#1#2}}}
