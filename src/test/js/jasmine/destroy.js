describe("destroy", function(){
    "use strict";
    beforeAll(function(done){
        Clotho.login("read", "read")
        .finally(function(){
            done();
        });
    });
    
    afterAll(function(done){
        Clotho.logout()
        .finally(function(){
            done();
        });
    });

    var ids;

    beforeEach(function(done){
        Clotho.create([
            {name:"Delete Me 1"},
            {name:"Delete Me 2"}
        ])
        .then(function(result){
            ids = result;
        })
        .finally(function(){
            done();
        });
    });

    afterEach(function(done){
        Clotho.destroy(ids)
        .finally(function(){
            done();
        });
    });

    it("destroys objects", function(done){
        Clotho.get(ids[0])
        .then(function(result){
            expect(result).not.toBeUndefined();
            return Clotho.delete(ids[0])
        })
        .then(function(){
            return Clotho.get(ids[0]);
        })
        .then(function(result){
            expect(result).toBeUndefined();
        })
        .finally(function(){
            done();
        });
    });

    describe("destroyAll", function(){
        it("destroys multiple objects", function(done){
            Clotho.get(ids)
            .then(function(result){
                expect(result[0]).not.toBeNull();
                expect(result[1]).not.toBeNull();
                return Clotho.delete(ids);
            })
            .then(function(){
               return Clotho.get(ids);
            })
            .then(function(result){
                expect(result[0]).toBeNull();
                expect(rseult[1]).toBeNull();
            })
            .finally(function(){
                done();
            });
        });
    });
});

