describe("create", function() {
    "use strict";

    beforeAll(function(done){
        Clotho.login("read", "read")
        .finally(function(){
            done();
        })
    })

    afterAll(function(done){
        Clotho.logout()
        .finally(function(){
            done();
        })
    })

    beforeEach(function(){
        toDelete = [];
    });

    afterEach(function(done){
        Clotho.destroy(toDelete)
        .finally(function(result){
            done();
        });
    });

    var toDelete;

    it("returns an id", function(done){
        Clotho.create({})
        .then(function(result){
            toDelete.push(result);
            expect(result).toEqual(jasmine.any(String));
        })
        .finally(function(){
            done();
        });
    });

    it("creates specified object in database", function(done){
        Clotho.create({"name":"Created Object"})
        .then(function(result){
            toDelete.push(result);
            return Clotho.get(result);
        })
        .then(function(result){
            expect(result.name).toEqual("Created Object");
        })
        .finally(function(){
            done();
        });
    });

    it("fails if id is preexisting", function(done){
        var id;
        Clotho.create({})
        .then(function(result){
            toDelete.push(result);
            id = result;
            return Clotho.create({"id":id, "name":"Second Attempt"});
        })
        .then(function(result){
            return Clotho.get(id);
        })
        .then(function(result){
            expect(result.name).toBeUndefined();
        })
        .finally(function(){
            done();
        });
    });

    it("returns undefined if it fails", function(done){
        var id;
        Clotho.create({})
        .then(function(result){
            toDelete.push(result);
            id = result;
            return Clotho.create({"id":id, "name":"Second Attempt"});
        })
        .then(function(result){
            expect(result).toBeUndefined();
        })
        .finally(function(){
            done();
        });
    });

    describe("createAll", function(){
        var ids;
        var data;
        beforeEach(function(done){
            var id;
            Clotho.create({})
            .then(function(result){
                id = result;
                return Clotho.create([
                    {"name":"first object"},
                    {"name":"second object"},
                    {"id":id, "name":"Second Attempt"}
                ]);
            })
            .then(function(result){
                ids = result;
                return Clotho.get(ids);
            })
            .then(function(result){
                data = result;
            })
            .finally(function(){
                done();
            });
        });

        it("can create multiple objects in db", function(){
            expect(data[0]).toBeTruthy();
            expect(data[1]).toBeTruthy();
        });

        it("returns ids in order", function(){
            expect(data[0].name).toEqual("first object");
            expect(data[1].name).toEqual("second object");
        });

        it("fails if id preexisting", function(){
            expect(data[2]).toBeNull();
        });

        it("returns null for failed item", function(){
            expect(ids[2]).toBeNull();
        });
    });
});

