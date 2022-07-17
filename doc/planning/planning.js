var tasks = [
{"startDate":new Date(2014,01,19),"endDate":new Date(2014,02,28),"taskName":"Rédaction article","status":"RUNNING"},
{"startDate":new Date(2014,02,28),"endDate":new Date(2014,02,28,12),"taskName":"Soumission article","status":"EVENT"},
{"startDate":new Date(2014,02,24),"endDate":new Date(2014,02,31),"taskName":"Plan manuscrit","status":"WAITING"},
{"startDate":new Date(2014,03,01),"endDate":new Date(2014,06,01),"taskName":"Rédaction manuscrit","status":"WAITING"},
{"startDate":new Date(2014,06,01),"endDate":new Date(2014,06,01,12),"taskName":"Soumission manuscrit","status":"EVENT"},
{"startDate":new Date(2014,02,24),"endDate":new Date(2014,04,30),"taskName":"Runs","status":"WAITING"},
];

var taskStatus = {
    "WAITING" : "bar",
    "RUNNING" : "bar-running",
    "DONE" : "bar-done", 
    "EVENT" : "events"
};

var taskNames = [];
for (var i = 0; i < tasks.length; i++)
  taskNames[i] = tasks[i].taskName;

var format = "%d/%m";

var gantt = d3.gantt().taskTypes(taskNames).taskStatus(taskStatus).tickFormat(format);

gantt(tasks);
