var canvas = document.getElementById("myCanvas");
var canvasWidth = canvas.width;
var canvasHeight = canvas.height;
var ctx = canvas.getContext("2d");
var canvasData = ctx.getImageData(0, 0, canvasWidth, canvasHeight);
var datefield = document.getElementById('date');
var fps = 50;
var offset = [0, 0];
var newOffset = [0, 0];
var focus = 'Earth';
var iter = 0;
var zoom = 50000;
var sidereal = 365.256366;
var au = 149597870700;
var cy2sec = 365.25 * 100 * 24 * 60 * 60 * 60;
var J2000 = Date.UTC(2000, 0, 1, 12, 0, 0);
var tMillisFromJ2000 = Date.now() - J2000;
var incMillis = fps * 5;
//1000 * 60 * 10;
var cMillis = 1000 * 60 * 60 * 24 * 365.25 * 100;
var lineWiggle = 0.25;
var lineWigglePixelX = lineWiggle * canvasWidth;
var lineWigglePixelY = lineWiggle * canvasHeight;
var tCenturiesFromJ2000 = 0;
var newZoom = 0;
var deltaMove = [0, 0];
var isMoving = false;
var G = 6.67408e-11;
var orbitDegreeMultiplier = 1;
var maxZoom = 50000;
var minZoom = 0.01;

window.bodies = [];
window.bodiesIndex = {};

function au2pixel(pos) {
    return Math.floor(pos * 100 * zoom);
}

function drawPixel(x, y, r, g, b, a) {

    plotPoint(x, y, [r, g, b, a]);
}

function drawLine(x0, y0, x1, y1, r, g, b, a) {
    color = [r, g, b, a];

    if (Math.abs(y1 - y0) < Math.abs(x1 - x0)) {
        if (x0 > x1)
            plotLineLow(x1, y1, x0, y0, color);
        else
            plotLineLow(x0, y0, x1, y1, color)

    } else {
        if (y0 > y1)
            plotLineHigh(x1, y1, x0, y0, color)
        else
            plotLineHigh(x0, y0, x1, y1, color)

    }
}

function plotPoint(x, y, color) {
    x = x + 200 - offset[0];
    y = y + 200 - offset[1];
    if (x < 0 || y < 0 || x > canvasWidth - 1 || y > canvasHeight - 1) {
        return;
    }
    var index = (x + y * canvasWidth) * 4;

    canvasData.data[index + 0] = color[0];
    canvasData.data[index + 1] = color[1];
    canvasData.data[index + 2] = color[2];
    canvasData.data[index + 3] = color[3];
}

function plotLineLow(x0, y0, x1, y1, color) {
    dx = x1 - x0;
    dy = y1 - y0;
    yi = 1;
    if (dy < 0) {
        yi = -1;
        dy = -dy;
    }

    D = 2 * dy - dx;
    y = y0;
    for (var x = x0; x <= x1; x++) {
        plotPoint(x, y, color);
        if (D > 0) {
            y = y + yi;
            D = D - 2 * dx;
        }
        D = D + 2 * dy;
    }
}

function plotLineHigh(x0, y0, x1, y1, color) {
    dx = x1 - x0;
    dy = y1 - y0;
    xi = 1;
    if (dx < 0) {
        xi = -1;
        dx = -dx;
    }

    D = 2 * dx - dy;
    x = x0;

    for (var y = y0; y <= y1; y++) {
        plotPoint(x, y, color);
        if (D > 0) {
            x = x + xi
            D = D - 2 * dy
        }
        D = D + 2 * dx
    }
}

function IEEERemainder(x, y) {
    var regularMod = x % y;
    if (isNaN(regularMod)) {
        return NaN;
    }
    if (regularMod == 0) {
        if (x < 0) {
            return -0;
        }
    }
    var alternativeResult;
    alternativeResult = regularMod - (Math.abs(y) * Math.sign(x));
    if (Math.abs(alternativeResult) == Math.abs(regularMod)) {
        var divisionResult = x / y;
        var roundedResult = Math.round(divisionResult);
        if (Math.abs(roundedResult) > Math.abs(divisionResult)) {
            return alternativeResult;
        } else {
            return regularMod;
        }
    }
    if (Math.abs(alternativeResult) < Math.abs(regularMod)) {
        return alternativeResult;
    } else {
        return regularMod;
    }
}

function timeConverter(UNIX_timestamp) {
    var a = new Date(UNIX_timestamp);
    var months = ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'];
    var year = a.getFullYear();
    var month = months[a.getMonth()];
    var date = a.getDate();
    var hour = a.getHours();
    var min = a.getMinutes();
    var sec = a.getSeconds();
    return date + ' ' + month + ' ' + year + ' ' + hour + ':' + min + ':' + sec;
}

function deg2rad(degrees) {
    return degrees * (Math.PI / 180);
}

function rad2deg(radians) {
    return radians / (Math.PI / 180);
}

function revsDayCount2revsCyDegrees(revsDay) {
    return revsDay * 360 * sidereal * 100
}

function updateCanvas() {
    ctx.putImageData(canvasData, 0, 0);
}

var bodiesConfig = [
    {
        name: 'Sol',
        mass: 1.9885e+30,
        pos: [0, 0, 0],
        color: [255, 255, 255, 255],
        parentBody: null,
        radius: 695508000

    },
    {
        name: 'Mercury',
        semiMajorAxis: 0.38709843,
        eccentricity: 0.20563661,
        inclination: 7.00559432,
        longMean: 252.25166724,
        longPeri: 77.45771895,
        longAscNode: 48.33961819,
        semiMajorAxisCy: 0.00000000,
        eccentricityCy: 0.00002123,
        inclinationCy: -0.00590158,
        longMeanCy: 149472.67486623,
        longPeriCy: 0.15940013,
        longAscNodeCy: -0.12214182,
        mass: 3.3011e23,
        color: [255, 255, 155, 255],
        parentBody: 'Sol',
        radius: 2439700
    },
    {
        name: 'Venus',
        semiMajorAxis: 0.72332102,
        eccentricity: 0.00676399,
        inclination: 3.39777545,
        longMean: 181.97970850,
        longPeri: 131.76755713,
        longAscNode: 76.67261496,
        semiMajorAxisCy: -0.00000026,
        eccentricityCy: -0.00005107,
        inclinationCy: 0.00043494,
        longMeanCy: 58517.81560260,
        longPeriCy: 0.05679648,
        longAscNodeCy: -0.27274174,
        mass: 4.8675e24,
        color: [255, 0, 0, 255],
        parentBody: 'Sol',
        radius: 6051800
    },
    {
        name: 'Earth',
        semiMajorAxis: 1.00000018,
        eccentricity: 0.01673163,
        inclination: -0.00054346,
        longMean: 100.46691572,
        longPeri: 102.93005885,
        longAscNode: -5.11260389,
        semiMajorAxisCy: -0.00000003,
        eccentricityCy: -0.00003661,
        inclinationCy: -0.01337178,
        longMeanCy: 35999.37306329,
        longPeriCy: 0.31795260,
        longAscNodeCy: -0.24123856,
        mass: 5.97237e24,
        color: [155, 155, 255, 255],
        parentBody: 'Sol',
        radius: 6378100
    },
    {
        name: 'Moon',
        semiMajorAxis: 0.00256955529,
        eccentricity: 0.0549,
        inclination: 5.145,
        longMean: 0,
        longPeri: 77.45771895,
        longAscNode: 48.33961819,
        semiMajorAxisCy: 0.00000000,
        eccentricityCy: 0.00002123,
        inclinationCy: -0.00590158,
        longMeanCy: 35999.37306329 * 12,
        longPeriCy: 0.15940013,
        longAscNodeCy: -0.12214182,
        mass: 7.342e22,
        color: [0, 255, 255, 255],
        parentBody: 'Earth',
        pos: [0, 0, 0],
        radius: 1738100
    },
    {
        name: 'Mars',
        semiMajorAxis: 1.52371243,
        eccentricity: 0.09336511,
        inclination: 1.85181869,
        longMean: -4.56813164,
        longPeri: -23.91744784,
        longAscNode: 49.71320984,
        semiMajorAxisCy: 0.00000097,
        eccentricityCy: 0.00009149,
        inclinationCy: -0.00724757,
        longMeanCy: 19140.29934243,
        longPeriCy: 0.45223625,
        longAscNodeCy: -0.26852431,
        parentBody: 'Sol',
        radius: 3396200
    },
    {
        name: 'Jupiter',
        semiMajorAxis: 5.20248019,
        eccentricity: 0.04853590,
        inclination: 1.29861416,
        longMean: 34.33479152,
        longPeri: 14.27495244,
        longAscNode: 100.29282654,
        semiMajorAxisCy: -0.00002864,
        eccentricityCy: 0.00018026,
        inclinationCy: -0.00322699,
        longMeanCy: 3034.90371757,
        longPeriCy: 0.18199196,
        longAscNodeCy: 0.13024619,
        parentBody: 'Sol',
        radius: 71492000
    },
    {
        name: 'Saturn',
        semiMajorAxis: 9.54149883,
        eccentricity: 0.05550825,
        inclination: 2.49424102,
        longMean: 50.07571329,
        longPeri: 92.86136063,
        longAscNode: 113.63998702,
        semiMajorAxisCy: -0.00003065,
        eccentricityCy: -0.00032044,
        inclinationCy: 0.00451969,
        longMeanCy: 1222.11494724,
        longPeriCy: 0.54179478,
        longAscNodeCy: -0.25015002,
        parentBody: 'Sol',
        radius: 60268000
    },
    {
        name: 'Uranus',
        semiMajorAxis: 19.18797948,
        eccentricity: 0.04685740,
        inclination: 0.77298127,
        longMean: 314.20276625,
        longPeri: 172.43404441,
        longAscNode: 73.96250215,
        semiMajorAxisCy: -0.00020455,
        eccentricityCy: -0.00001550,
        inclinationCy: -0.00180155,
        longMeanCy: 428.49512595,
        longPeriCy: 0.09266985,
        longAscNodeCy: 0.05739699,
        parentBody: 'Sol',
        radius: 25559000
    },
    {
        name: 'Neptune',
        semiMajorAxis: 30.06952752,
        eccentricity: 0.00895439,
        inclination: 1.77005520,
        longMean: 304.22289287,
        longPeri: 46.68158724,
        longAscNode: 131.78635853,
        semiMajorAxisCy: 0.00006447,
        eccentricityCy: 0.00000818,
        inclinationCy: 0.00022400,
        longMeanCy: 218.46515314,
        longPeriCy: 0.01009938,
        longAscNodeCy: -0.00606302,
        parentBody: 'Sol',
        radius: 24764000
    },
    {
        name: 'Pluto',
        semiMajorAxis: 39.48686035,
        eccentricity: 0.24885238,
        inclination: 17.14104260,
        longMean: 238.96535011,
        longPeri: 224.09702598,
        longAscNode: 110.30167986,
        semiMajorAxisCy: 0.00449751,
        eccentricityCy: 0.00006016,
        inclinationCy: 0.00000501,
        longMeanCy: 145.18042903,
        longPeriCy: -0.00968827,
        longAscNodeCy: -0.00809981,
        parentBody: 'Sol',
        radius: 1195000
    },
    {
        name: 'ISS',
        semiMajorAxis: (403156.5655 + 6378100) / au,
        eccentricity: 0.0003885,
        inclination: 51.6373,
        longMean: 153.1203, // degrees past starting?
        longPeri: 206.9748,
        longAscNode: 238.6885,
        longMeanCy: revsDayCount2revsCyDegrees(15.53729445),
        mass: 10000,
        color: [255, 255, 155, 255],
        parentBody: 'Earth',
        radius: 50
    },
    {
        name: 'SS 1',
        semiMajorAxis: 6378100 / au,
        eccentricity: 0,
        inclination: 0.6373,
        longMean: 0, // degrees past starting?
        longPeri: 206.9748,
        longAscNode: 238.6885,
        longMeanCy: revsDayCount2revsCyDegrees(1),
        mass: 10000,
        color: [255, 255, 155, 255],
        parentBody: 'Earth',
        radius: 50,
        //history: []
    }
];

for (var i = 0; i < bodiesConfig.length; i++) {
    window.bodiesIndex[bodiesConfig[i].name] = i;
}

//http://www.jgiesen.de/astro/astroJS/siderealClock/sidClock.js
function GM_Sidereal_Time(jd) {
    var t_eph, ut, MJD0, MJD;

    MJD = jd - 2400000.5;
    MJD0 = Math.floor(MJD);
    ut = (MJD - MJD0) * 24.0;
    t_eph = (MJD0 - 51544.5) / 36525.0;
    return 6.697374558 + 1.0027379093 * ut + (8640184.812866 + (0.093104 - 0.0000062 * t_eph) * t_eph) * t_eph / 3600.0;
}

function LM_Sidereal_Time(jd, longitude) {
    var GMST = GM_Sidereal_Time(jd);
    var LMST = 24.0 * frac((GMST + longitude / 15.0) / 24.0);
    return (LMST);
}

function frac(X) {
    X = X - Math.floor(X);
    if (X < 0) X = X + 1.0;
    return X;
}

function HoursMinutesSeconds(time) {

    var h = Math.floor(time);
    var min = Math.floor(60.0 * frac(time));
    var secs = Math.round(60.0 * (60.0 * frac(time) - min));


    return 15 * (h + min / 60 + secs / 3600)

}

function JulDay(d, m, y, u) {
    if (y < 1900) y = y + 1900
    if (m <= 2) {
        m = m + 12;
        y = y - 1
    }
    A = Math.floor(y / 100);
    JD = Math.floor(365.25 * (y + 4716)) + Math.floor(30.6001 * (m + 1)) + d - 13 - 1524.5 + u / 24.0;
    return JD
}

var bodiesMethods = {
    siderealTime: function (longitude, time) {

        if (time === undefined) {
            var time = new Date(this.time + J2000);
            time = new Date(time.valueOf() + time.getTimezoneOffset() * 60000);
        }

        UT = time.getHours() + time.getMinutes() / 60 + time.getSeconds() / 3600;

        //console.log(p);
        JD = JulDay(time.getDate(), time.getMonth() + 1, time.getYear(), UT);
        return LM_Sidereal_Time(JD, longitude) / 24 * 360;
    },
    rApoapsis: function () {
        return this.semiMajorAxis * (1 + this.eccentricity);
    },
    rPeriapsis: function () {
        return this.semiMajorAxis * (1 - this.eccentricity);
    },
    rCurrent: function () {
        var x = this.relativePos[0];
        var y = this.relativePos[1];
        var z = this.relativePos[2];

        return Math.sqrt(x * x + y * y + z * z) * au;
    },
    trueAnomaly: function () {
        return IEEERemainder(IEEERemainder(this.longMeanCy * this.timeCy, 360) + this.longMean, 360);
    },
    eccentricAnomaly: function () {
        var radians = Math.PI / 180 * this.trueAnomaly();
        return Math.acos((this.eccentricity + Math.cos(radians)) / (1 + this.eccentricity * Math.cos(radians)));
    },
    flightPathAngle: function () {
        var a = this.velocity();
        return Math.atan(this.eccentricity * Math.sin(a) / (1 + this.eccentricity * Math.cos(a)));
    },
    flightPathAngle2: function () {
        var v = this.velocity();
        var r = this.rCurrent();

        return Math.asin((v * r) / (Math.abs(v) * Math.abs(r)));
    },
    velocity: function () {
        var a = this.semiMajorAxis;
        return Math.sqrt((this.mass + this.parent().mass) * G * (2 / this.rCurrent() - 1 / (a * au)));
    },
    // this is the next 360 degrees of path
    compute: function () {
        if (this.parentBody === null)
            return;
        var out = [];
        for (var i = 0; i < 360 * orbitDegreeMultiplier; i++) {
            var L = i - 180;
            this.keplar(L / orbitDegreeMultiplier);
            out.push(this.pos);
        }
        this.computed = out;
    },
    keplar: function (L) {
        if (this.parentBody === null) {
            return;
        }
        var a = this.semiMajorAxis // + body.semiMajorAxisCy * tCenturiesFromJ2000;
        var e = this.eccentricity // + body.eccentricityCy * tCenturiesFromJ2000;
        var i = this.inclination // + body.inclinationCy * tCenturiesFromJ2000;
        if (L === undefined)
            L = IEEERemainder(IEEERemainder(this.longMeanCy * this.timeCy, 360) + this.longMean, 360);
        var p = this.longPeri // + body.longPeriCy * tCenturiesFromJ2000;
        var W = this.longAscNode // + body.longAscNodeCy * tCenturiesFromJ2000;

        var M = L - p; // p is the longitude of periapsis
        var w = p - W; // W is the longitude of the ascending node

        w = deg2rad(w);
        M = deg2rad(M);
        i = deg2rad(i);

        var E = M;
        var iters = 0;
        while (true) {
            iters++;
            var dE = (E - e * Math.sin(E) - M) / (1 - e * Math.cos(E));
            E -= dE;
            if (Math.abs(dE) < 1e-6 || iters > 500) break;
        }

        var P = a * (Math.cos(E) - e);

        var Q = a * Math.sin(E) * Math.sqrt(1 - Math.pow(e, 2));

        // rotate by argument of periapsis
        var x = Math.cos(w) * P - Math.sin(w) * Q;
        var y = Math.sin(w) * P + Math.cos(w) * Q;
        // rotate by inclination
        var z = Math.sin(i) * x;
        x = Math.cos(i) * x;
        // rotate by longitude of ascending node
        var xtemp = x;
        x = Math.cos(W) * xtemp - Math.sin(W) * y;
        y = Math.sin(W) * xtemp + Math.cos(W) * y;

        if (this.history !== undefined && this.relativePos !== undefined) {
            this.history.push(this.relativePos.slice(0));
            if (this.history.length > 500) {
                this.history.shift();
            }
        }
        this.pos = [x, y, z];
        this.relativePos = [x, y, z];
    },
    calculateAbsolutePosition: function () {
        if (this.parentBody === null) {
            return;
        }
        this.pos[0] += this.parent().pos[0];
        this.pos[1] += this.parent().pos[1];
        this.pos[2] += this.parent().pos[2];
    },
    periApoFromRVZ: function (radius, velocity, zenithAngle, body2) {
        M = this.mass;
        if (body2 !== undefined)
            M += body2.mass;

        C = 2 * G * M / (radius * Math.pow(velocity, 2));

        oneMinC = 1 - C;
        cSquared = C * C;

        rPA1 = (-C + Math.sqrt(cSquared - 4 * oneMinC * -Math.pow(Math.sin(deg2rad(zenithAngle)), 2))) / (2 * oneMinC);
        rPA2 = (-C - Math.sqrt(cSquared - 4 * oneMinC * -Math.pow(Math.sin(deg2rad(zenithAngle)), 2))) / (2 * oneMinC);

        if (rPA1 < rPA2) {
            return {
                rPerigee: rPA1 * radius,
                rApogee: rPA2 * radius
            }
        }
        return {
            rApogee: rPA1 * radius,
            rPerigee: rPA2 * radius
        }
    },
    periApoAltitudeFromRVZ: function (altitude, velocity, zenithAngle, body2) {
        result = this.periApoFromRVZ(altitude + this.radius, velocity, zenithAngle, body2);
        result.rPerigee -= this.radius;
        result.rApogee -= this.radius;
        return result;
    },
    eccentricityFromRVZ: function (radius, velocity, zenithAngle, body2) {
        M = this.mass;
        if (body2 !== undefined)
            M += body2.mass;

        C = 2 * G * M / (radius * Math.pow(velocity, 2));

        sin2 = Math.pow(Math.sin(deg2rad(zenithAngle)), 2);
        cos2 = Math.pow(Math.cos(deg2rad(zenithAngle)), 2);
        return Math.sqrt(Math.pow((radius * (velocity * velocity) / (G * M) - 1), 2) * sin2 + cos2);
    },
    anglePeriFromRVZ: function (radius, velocity, zenithAngle, body2) {
        M = this.mass;
        if (body2 !== undefined)
            M += body2.mass;

        C = 2 * G * M / (radius * Math.pow(velocity, 2));
        vSquared = velocity * velocity;
        zen = deg2rad(zenithAngle);
        GM = G * M;
        sin2 = Math.pow(Math.sin(zen), 2);
//(6,628,140 × 7,9002 / 3.986005×1014) × sin(89) × cos(89)
// 	       / [(6,628,140 × 7,9002 / 3.986005×1014) × sin2(89) - 1]

        tan = (radius * vSquared / GM) * Math.sin(zen) * Math.cos(zen) / ((radius * vSquared / GM) * sin2 - 1);

        return rad2deg(Math.atan(tan));

    },
    semiMajorAxisFromRV: function (radius, velocity, body2) {
        M = this.mass;
        if (body2 !== undefined)
            M += body2.mass;

        C = 2 * G * M / (radius * Math.pow(velocity, 2));
        vSquared = velocity * velocity;
        GM = G * M;

        return 1 / (2 / radius - vSquared / GM)
    },
    inclinationFromLatAzimuth: function (latitude, azimuth) {
        return rad2deg(Math.acos(Math.cos(deg2rad(latitude)) * Math.sin(deg2rad(azimuth))));
    },
    angularDistanceOrbital: function (latitude, azimuth) {
        return rad2deg(Math.atan(Math.tan(deg2rad(latitude)) / Math.cos(deg2rad(azimuth))));
    },
    angularDistanceEquatorial: function (latitude, azimuth) {
        return rad2deg(Math.atan(Math.sin(deg2rad(latitude)) * Math.tan(deg2rad(azimuth))))
    },
    longitudeAscNodeFrom: function (angularDistanceEquitorial, insertionLongitude) {
        return this.siderealTime(insertionLongitude - angularDistanceEquitorial);

    },
    argPerigee: function (anglePeri, angularDistance) {
        // angularDistanceOrbital
        return angularDistance - anglePeri;
    },
    parent: function () {
        if (window.bodiesIndex[this.parentBody] === this) return null;
        return window.bodies[window.bodiesIndex[this.parentBody]];
    },
    drawOrbit: function () {
        if (this.computed === undefined)
            return;
        for (var i = 0; i < this.computed.length; i++) {
            last = i - 1;
            if (last < 0)
                last += this.computed.length;

            x0 = au2pixel(this.computed[i][0] + this.parent().pos[0]);
            y0 = au2pixel(this.computed[i][1] + this.parent().pos[1]);

            x1 = au2pixel(this.computed[last][0] + this.parent().pos[0]);
            y1 = au2pixel(this.computed[last][1] + this.parent().pos[1]);

            if (x0 + lineWigglePixelX + 200 - offset[0] < 0 && x1 + lineWigglePixelX + 200 - offset[0] < 0)
                continue;
            if (y0 + lineWigglePixelY + 200 - offset[1] < 0 && y1 + lineWigglePixelY + 200 - offset[1] < 0)
                continue;
            if (x0 - lineWigglePixelX + 200 - offset[0] > 600 && x1 - lineWigglePixelX + 200 - offset[0] > 600)
                continue;
            if (y0 - lineWigglePixelY + 200 - offset[1] > 600 && y1 + lineWigglePixelY + 200 - offset[1] > 600)
                continue;
            drawLine(x0, y0, x1, y1, 155, 155, 155, 255)
        }
        if (this.history !== undefined)
            for (var i = 0; i+1 < this.history.length; i++) {


                x0 = au2pixel(this.history[i][0] + this.parent().pos[0]);
                y0 = au2pixel(this.history[i][1] + this.parent().pos[1]);

                x1 = au2pixel(this.history[i+1][0] + this.parent().pos[0]);
                y1 = au2pixel(this.history[i+1][1] + this.parent().pos[1]);

                if (x0 + lineWigglePixelX + 200 - offset[0] < 0 && x1 + lineWigglePixelX + 200 - offset[0] < 0)
                    continue;
                if (y0 + lineWigglePixelY + 200 - offset[1] < 0 && y1 + lineWigglePixelY + 200 - offset[1] < 0)
                    continue;
                if (x0 - lineWigglePixelX + 200 - offset[0] > 600 && x1 - lineWigglePixelX + 200 - offset[0] > 600)
                    continue;
                if (y0 - lineWigglePixelY + 200 - offset[1] > 600 && y1 + lineWigglePixelY + 200 - offset[1] > 600)
                    continue;
                drawLine(x0, y0, x1, y1, 200, 155, 200, 255)
            }
    },
    drawBody: function (sphere) {
        var color = [255, 255, 255, 255];
        if (this.hasOwnProperty('color')) {
            color = this.color;
        }
        if (sphere) {
            x = au2pixel(this.pos[0]) + 200 - offset[0];
            y = au2pixel(this.pos[1]) + 200 - offset[1];
            if (x > 0 && x < 600 && y > 0 && y < 600) {
                ctx.beginPath();

                ctx.fillStyle="#"+color[0].toString(16)+color[1].toString(16)+color[2].toString(16)+"66";
                ctx.arc(x, y, au2pixel(this.radius / au),0,2*Math.PI);
                ctx.fill();
            }
        } else {

            plotPoint(au2pixel(this.pos[0]), au2pixel(this.pos[1]), color);
            plotPoint(au2pixel(this.pos[0]) + 1, au2pixel(this.pos[1]), color);
            plotPoint(au2pixel(this.pos[0]), au2pixel(this.pos[1]) + 1, color);
            plotPoint(au2pixel(this.pos[0]) - 1, au2pixel(this.pos[1]), color);
            plotPoint(au2pixel(this.pos[0]), au2pixel(this.pos[1]) - 1, color);
        }
    }

};

// init the bodies
for (var i = 0; i < bodiesConfig.length; i++) {
    window.bodies[i] = {};
    for (var j = 0; j < Object.keys(bodiesMethods).length; j++) {
        var k = Object.keys(bodiesMethods)[j];
        window.bodies[i][k] = bodiesMethods[k];
    }
    for (var j = 0; j < Object.keys(bodiesConfig[i]).length; j++) {
        var k = Object.keys(bodiesConfig[i])[j];
        window.bodies[i][k] = bodiesConfig[i][k];
    }

}

for (var i = 0; i < window.bodies.length; i++) {
    window.bodies[i].timeCy = (tMillisFromJ2000 + (incMillis * iter++)) / cMillis;
    window.bodies[i].time = 0;
    window.bodies[i].compute();
}

var targetBodyStats = 12;
var tMission = 0;

function tMission2velocity() {
    i  = tMission/1000/600;
    if (i > 1) i = 1
    return i * 7500;
}

function earthAtmosphericDensity(km) {
    var density = [1.17E+00,9.48E-02,4.07E-03,3.31E-04,1.69E-05,5.77E-07,1.7E-08,2.96E-09,9.65E-10,3.9E-10,1.75E-10,8.47E-11,4.31E-11,2.3E-11,1.27E-11,7.22E-12,4.21E-12,2.5E-12,1.51E-12,9.2E-13,5.68E-13,3.54E-13,2.23E-13,1.42E-13,9.2E-14,6.03E-14,4.03E-14,2.75E-14,1.93E-14,1.39E-14,1.03E-14,7.9E-15,6.24E-15,5.06E-15,4.21E-15,3.58E-15,3.09E-15,2.7E-15,2.39E-15,2.13E-15,1.91E-15,1.73E-15,1.56E-15,1.42E-15,1.3E-15,1.18E-15];
    item = Math.floor(km/20);
    if (item < 0) { item = 0 };
    if (item >= density.length) { return 0 };
    return density[item];
}


function drag(dragCoe, surfaceArea, radius, velocity) {
    return (Math.PI * dragCoe * surfaceArea * earthAtmosphericDensity(radius/1000) * radius * velocity) / 1000;
}
var text_ = 2;
function text()
{
    text_ += 12;
    return text_;
}
var totalDelta = 0;

var ss = window.bodies[12];
var radius = ss.rCurrent();
var updateFn = function () {
    text_ = 2;
    ctx.beginPath();
    ctx.rect(0, 0, 600, 600);
    ctx.fillStyle = "black";
    ctx.fill();

    ctx.font = "12px Arial";
    ctx.fillStyle = '#DDDDFF';
    canvasData = ctx.getImageData(0, 0, canvasWidth, canvasHeight);
    tCenturiesFromJ2000 = (tMillisFromJ2000 + (incMillis * iter++)) / cMillis;
    datefield.innerHTML = timeConverter(J2000 + cMillis * tCenturiesFromJ2000);
    tMission = (incMillis * iter++);

    for (var b = 0; b < window.bodies.length; b++) {
        window.bodies[b].drawOrbit();
        window.bodies[b].timeCy = tCenturiesFromJ2000;
        window.bodies[b].time = cMillis * tCenturiesFromJ2000;
        window.bodies[b].keplar();
        window.bodies[b].calculateAbsolutePosition();

        if (focus === window.bodies[b].name) {
            offset = [au2pixel(window.bodies[b].pos[0]), au2pixel(window.bodies[b].pos[1])];
        }

        window.bodies[b].drawBody(false);
    }

    updateCanvas();
    for (var b = 0; b < window.bodies.length; b++) {
        window.bodies[b].drawBody(true);
    }
    ctx.fillStyle="#FFFFFF";
    if (newZoom > 0) {
        zoom = newZoom;
        offset = [newOffset[0],newOffset[1]];

        newZoom = null;
    }
    var earth = window.bodies[3];
    ctx.fillText(window.bodies[targetBodyStats].name + " FPA: " + window.bodies[targetBodyStats].flightPathAngle().toFixed(6) + 'deg', 10, text());
    ctx.fillText(window.bodies[targetBodyStats].name + " TA: " + window.bodies[targetBodyStats].trueAnomaly().toFixed(2) + 'deg', 10, text());

    velocity = tMission2velocity();
    var i = tMission/1000/600;
    if (i < 1) {
        totalDelta += velocity * incMillis/1000*0.5;

    } else {
        period = (24 * 60 * 60) / (Math.sqrt (Math.pow(ss.semiMajorAxis * au,3) / (G * earth.mass)) * 2 * Math.PI);
        ss.longMeanCy = revsDayCount2revsCyDegrees(period);
    }

    zenith = 89;
    radiusD = radius + totalDelta;
    //radiusD = radius + 400000;
    //velocity = 6000;
    parameters = earth.periApoFromRVZ(radiusD, velocity, zenith);


    //fGrav = ss.mass * 9.80665;

    ss.semiMajorAxis = earth.semiMajorAxisFromRV(radiusD, velocity) / au;
    ss.eccentricity = earth.eccentricityFromRVZ(radiusD, velocity, zenith);




    ss.meanLong = 180;

    ss.longPeri = earth.anglePeriFromRVZ(radiusD, velocity, zenith);
    tEcc = 0;//earth.longitudeAscNodeFrom(earth.angularDistanceEquatorial(0, 90), 0);

    ss.longAscNode = tEcc;


    GM = earth.mass * G;

    test  = parameters.rApogee - earth.radius - totalDelta;


    ctx.fillText("apogee: " + ((parameters.rApogee - earth.radius) / 1000).toFixed(2) + 'km', 10, text());
    ctx.fillText("perigee: " + ((parameters.rPerigee - earth.radius) / 1000).toFixed(2) + 'km', 10, text());

    ctx.fillText("altitude: " + ((ss.rCurrent() - earth.radius) / 1000).toFixed(2) + 'km', 10, text());
    ctx.fillText("velocity: " + ss.velocity().toFixed(2) + '', 10, text());
    ctx.fillText("change in radius from acceleration: " + totalDelta.toFixed() + '', 10, text());
    ctx.fillText("semiMajor: " + (ss.semiMajorAxis * au / 1000).toFixed(2) + 'km', 10, text());
    ctx.fillText("eccentricity: " + ss.eccentricity.toFixed(3), 10, text());
    ctx.fillText("z: " + test.toFixed(), 10, text());
    ctx.fillText("period: " + ss.longMeanCy.toFixed(6), 10, text());




    ctx.fillText("circle: " + au2pixel(earth.pos[0]) + '.' + au2pixel(earth.pos[1]) + '.' + au2pixel(earth.radius / au), 10, text());

    tt = ss.time;
    ss.compute();
    ss.time = tt;


    /* var p = window.bodies[targetBodyStats].siderealTime();//window.bodies[targetBodyStats].argPerigee(25.794, window.bodies[targetBodyStats].angularDistanceOrbital(32, 86));

     */
    setTimeout(updateFn, fps)
}
updateFn();

window.addEventListener('wheel', function (e) {
    e.preventDefault();
    if (newZoom === 0)
        newZoom = zoom;
    newOffset = offset;
    if (e.deltaY < 0 && zoom > minZoom) {
        newZoom += 0.03 * e.deltaY * newZoom;
    }
    if (e.deltaY > 0 && zoom < maxZoom) {
        newZoom += 0.03 * e.deltaY * newZoom;
    }
    newOffset = [Math.floor(offset[0] * newZoom/zoom),Math.floor(offset[1] * newZoom/zoom)];
});
window.addEventListener('mousedown', function (e) {
    deltaMove = [e.pageX + offset[0], e.pageY + offset[1]];
    isMoving = true;
});
window.addEventListener('mousemove', function (e) {
    if (isMoving) {
        offset = [deltaMove[0] - e.pageX, deltaMove[1] - e.pageY]
    }
});
window.addEventListener('mouseup', function (e) {

    isMoving = false;
});
