const readline = require('readline')

// while(true){
//     console.log('hello world')
// }
console.log('hello world, from mocksatjs')

let ω = []
let b = []

const rl = readline.createInterface({
    input: process.stdin,
    output: process.stdout
})
rl.on('line', input => {
    if (input.length < 4 || input.slice(0, 3) != ">>>") return
    input = input.slice(3)
    let key = input[0]
    let value = input.slice(1)
    if (key == 'ω') {
        ω = JSON.parse(value)
    } else if (key == 'b') {
        b = JSON.parse(value)
    } else if (key == '?') {
        let control = ω.map(v => -v * 0.001)
        control = JSON.stringify(control)
        console.log(`>>>M${control}`)
    }
    // console.log(ω, b)
    // console.log(key, value)
})

rl.on('close', e => {
    console.log(`eee${e}`)
})
