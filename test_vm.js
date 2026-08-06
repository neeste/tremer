const vm = require('vm');
const fs = require('fs');

const sandbox = {
    console: console,
    setTimeout: setTimeout,
    fetch: fetch,
    Module: {
        locateFile: function(path) { return './wasm/' + path; }
    },
    // mock DOM and browser environments if needed
    window: {},
    document: { currentScript: { src: './wasm/tremer.js' } }
};

sandbox.globalThis = sandbox;
vm.createContext(sandbox);

const wasmCode = fs.readFileSync('./wasm/tremer.js', 'utf8');
vm.runInContext(wasmCode, sandbox);

setTimeout(() => {
    console.log("Keys in Module:", Object.keys(sandbox.Module));
    console.log("Is FS in Module?", !!sandbox.Module.FS);
}, 500);
