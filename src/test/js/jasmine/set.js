describe("set", function(){
    "use strict";

    beforeAll(function(done){
        Clotho.login("read","read")
        .then(function(){
            done();
        });
    });

    afterAll(function(done){
        Clotho.logout()
        .then(function(){
            done();
        });
    });

    var ids;
    beforeEach(function(done){
        Clotho.create([{"field": "value"}, {}])
        .then(function(result){
            ids = result;
        })
        .finally(function(){
            done();
        });
    });

    afterEach(function(done){
        Clotho.destroy(ids)
        .finally(function(result){
            done();
        });
    });


    it("can add new fields", function(done){
        Clotho.set({"id":ids[1], "new field":"new value"})
        .then(function(result){
            return Clotho.get(ids[1]);
        })
        .then(function(result){
            expect(result["new field"]).toEqual("new value");
        })
        .finally(function(){
            done();
        });
    });

    it("can change preexisting fields", function(done){
        Clotho.set({"id":ids[0], "field":"different value"})
        .then(function(result){
            return Clotho.get(ids[0]);
        })
        .then(function(result){
            expect(result["field"]).toEqual("different value");
        })
        .finally(function(){
            done();
        });
    });
    
    it("creates nonexistent objects", function(done){
        Clotho.set({"id":"does not exist", "field":"value"})
        .then(function(result){
            return Clotho.get("does not exist");
        })
        .then(function(result){
            expect(result.id).toEqual("does not exist");
            expect(result.field).toEqual("value");
        })
        .finally(function(){
            done();
        });
    });

    describe("setAll", function(){
        it("can edit multiple objects", function(done){
            Clotho.set([
                {"id":ids[0], "field":"different value"}, 
                {"id":ids[1], "new field":"new value"}
            ])
            .then(function(result){
                return Clotho.get(ids);
            })
            .then(function(results){
                expect(results[0]["field"]).toEqual("different value");
                expect(results[1]["new field"]).toEqual("new value");
            })
            .finally(function(){
                done();
            });
        });

        it("creates nonexistent objects", function(done){
            Clotho.set([{id:"does not exist", "field":"value"}])
            .then(function(result){
                return Clotho.get(["does not exist"]);
            })
            .then(function(results){
                expect(results[0].id).toEqual("does not exist");
                expect(results[0].field).toEqual("value");
            })
            .finally(function(){
                done();
            });
        });
    });
});
