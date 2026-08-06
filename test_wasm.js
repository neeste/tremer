const fs = require('fs');
const Module = require('./wasm/tremer.js');

Module.onRuntimeInitialized = () => {
    try {
        console.log("Runtime initialized");
        // write a dummy file
        const data = fs.readFileSync('Neely_project/Neely_s10.txt');
        Module.FS.writeFile('Neely_s10.txt', data);
        console.log("Calling main...");
        Module.callMain(['Neely_s10.txt']);
        console.log("Main finished successfully");
        
        console.log(Module.FS.readdir('.'));
    } catch (e) {
        console.error("Error caught:", e);
    }
};
