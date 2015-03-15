describe("Clotho 3 API", function() {
    it("exists", function (){
        expect(Clotho).toBeDefined();
    });
});

describe("Q API", function(){
    it("exists", function (){
        expect(Q).toBeDefined();
    });
});

describe("Clotho Server", function(){
    it("exists", function (done){
        var errorHandler = function(label){
            return function(e){
                expect(label).toBeNull();
                expect(e).toBeNull();
                done();
            }
        }

        var socket = new WebSocket("wss://localhost:8443/websocket");
        socket.onopen = function (){
            expect(true).toBeTruthy();
            done();
        }

        socket.onerror = errorHandler('error');
        socket.onclose = errorHandler('close');
    });
});
