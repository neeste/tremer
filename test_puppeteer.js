const puppeteer = require('puppeteer');
const fs = require('fs');

(async () => {
    try {
        const browser = await puppeteer.launch({ args: ['--no-sandbox'] });
        const page = await browser.newPage();
        
        page.on('console', msg => console.log('BROWSER CONSOLE:', msg.text()));
        page.on('pageerror', error => console.log('BROWSER ERROR:', error.message));
        
        await page.goto('file://' + __dirname + '/admin/index.html', { waitUntil: 'networkidle0' });
        
        // upload file
        const fileInput = await page.$('#fileInput');
        await fileInput.uploadFile(__dirname + '/Neely_project/Neely_s10.txt');
        
        // click generate
        console.log("Clicking Generate...");
        await page.evaluate(() => processFile());
        
        // wait a bit
        await new Promise(r => setTimeout(r, 2000));
        
        await browser.close();
        console.log("Done");
    } catch (e) {
        console.error("Puppeteer error:", e);
    }
})();
