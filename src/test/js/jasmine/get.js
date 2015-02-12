describe("get", function(){
    "use strict";
    var failSpecOnError = function(error){
        expect(error).toBeNull();
    }
    describe("under standard operation", function(){
        beforeEach(function(done) {
            //Q and Jasmine this-scoping don't play together
            var thiz = this;
            Clotho.get("org.clothocad.model.Part")
            .then(function(data){
                thiz.data = data;
            })
            .catch(failSpecOnError)
            .finally(function(){ 
                done();
            });
        });

        it("retrieves object data", function(){
            expect(this.data.name).toEqual("Part");
        });

        it("provides permission info", function(){
            expect(this.data.$$permissions).toBeDefined();
        });
    });

    it("returns nothing if nonexistent", function(done) {
        Clotho.get("nonexistent id")
        .then(function(data) {
            expect(data).toBeUndefined();
        })
        .catch(failSpecOnError)
        .finally(function(){done();});
    });

    it("can be muted");
    it("reports failure if nonexistent");
});

describe("getAll", function(){
    var failSpecOnError = function(error){
        expect(error).toBeNull();
    }
    beforeEach(function(done) {
        var thiz = this;
        Clotho.get(["org.clothocad.model.Part", "invalid id"])
        .then(function(response) {
            thiz.data = response;
        })
        .catch(failSpecOnError)
        .finally(function(){done();});
    });

    it("retrieves object data for a list of ids, in order", function(){
        expect(this.data[0].name).toEqual("Part");
    });

    it("returns null for nonexistent items", function() {
        expect(this.data[1]).toBeNull();
    });

    it("provides permission info", function() {
        expect(this.data[0].$$permissions).toBeDefined();
    });
});


