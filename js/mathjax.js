MathJax = {
    tex: {
        packages: ['base'],
        inlineMath: [
            ['$', '$']
            ['\\(', '\\)']
        ],
        displayMath: [
            ['$$', '$$'],
            ['\\[', '\\]']
        ],
        processEscapes: true,
        processEnvironments: true,
        processRefs: true,
        digits: /^(?:[0-9]+(?:\{,\}[0-9]{3})*(?:\.[0-9]*)?|\.[0-9]+)/,
        tags: 'none',
        tagSide: 'right',
        tagIndent: '0.8em',
        useLabelIds: true,
        multlineWidth: '85%',
        maxMacros: 1000,
        maxBuffer: 5 * 1024,
    }
};

(() => {
    var head = document.getElementsByTagName("head")[0], script;
    script = document.createElement("script");
    script.type = "text/javascript";
    script.src  = "https://cdn.jsdelivr.net/npm/mathjax@3/es5/tex-mml-chtml.js";
    head.appendChild(script);

    let promise = Promise.resolve();
    let typeset = (code) => {
        promise = promise.then(() => {
            let modify = code();
            return MathJax.typesetPromise(modify);
        }).catch((err) => {
            console.log('Typeset failed: ' + err.message);
            MathJax.texReset();
            MathJax.typesetClear();
            return;
        });
        return promise;
    };

    let mod = false;
    let observer = new MutationObserver(function(mutations){
        mutations.forEach(function(mutation) {
            mod = true;
        });
    })

    observer.observe(document.getElementById("left"), {
        characterData: true,
        subtree: true
    });

    setInterval(() => {
        if (mod) {
            typeset(() => {
                const right = document.getElementById("right");
                const left = document.getElementById("left");
                right.innerHTML = left.innerHTML.slice(0);
                return [right];
            });

            mod = false;
        }
    }, 200);
})();
