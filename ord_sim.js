createChart();

async function createChart() {
    const data = await getData();
    var ctx = document.getElementById('chart').getContext('2d');
    var myChart = new Chart(ctx, {
        type: 'scatter',
        data: {
            labels: data.xs,
            datasets: [{
                label: 'Voltage Over Time',
                data: data.ys,
                backgroundColor: 
                    'rgba(153, 102, 255, 0.2)',
                borderColor: 
                    'rgba(153, 102, 255, 0.5)',
                borderWidth: 1
            }]
        },
        options: {
            scales: {
                y: {
                    beginAtZero: true
                }
            }
        }
    });
}

async function getData() {
    var xs = [];
    var ys = [];

    const response = await fetch("Docs/ts.csv")
    const data = await response.text();

    const table = data.split('\n').slice(1);
    table.forEach(row => {
        const columns = row.split(',');
        const time = columns[0];
        xs.push(time);
        const voltage = columns[1];
        ys.push(voltage);
        //console.log(time, voltage)
    });
    return { xs, ys };
}

// xs = [];
// ys = [];

// time = Math.floor(t-(ft-beatssave*CL));
// xs.push(time);
// volt = Math.floor(v);
// ys.push(volt);