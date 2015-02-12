//not logged in
//proper permission info in get & query objects
//no permissions
//read permissions
//create without permissions
//create multiple without permissions
//write permissions
//run permissions
//owner permissions
//set w/o permissions
//set multiple w/o permissions

describe("authorization", function(){
    var publicObjectId;
    var privateObjectId; 
    var publicFunctionId;
    var privateFunctionId;

    var publicObject = "publicObject";
    var privateObject = "privateObject";
    var publicFunction = "publicFunction";
    var privateFunction = "privateFunction";

    var objectCopies;

    function resolveId(string){
        switch (string){
            case "publicObject":
                return publicObjectId;
            case "privateObject":
                return privateObjectId;
            case "publicFunction":
                return publicFunctionId;
            case "privateFunction":
                return privateFunctionId;
        }
    }

    beforeAll(function(done){
        Clotho.login("owner", "owner")
        .then(function(){
            return Clotho.create([
                {name:"Private Object", schema:"org.clothocad.model.Institution"},
                {name:"Public Object", schema:"org.clothocad.model.Institution"},
                {name:"Private Function", schema:"org.clothocad.core.util.SecurityTester"},
                {name:"Public Function", schema:"org.clothocad.core.util.SecurityTester"},
            ]);
        })
        .then(function(result){
            publicObjectId = result[0];
            privateObjectId = result[1];
            publicFunctionId = result[2];
            privateFunctionId = result[3];
            var promises = [];

            promises.push(Clotho.grant("read", result, ["read"], []));
            promises.push(Clotho.grant("write", result, ["write"], []));
            promises.push(Clotho.grant("run", result, ["run"], []));
            promises.push(Clotho.grant(null, [publicObjectId, publicFunctionId], ["public"], []));

            return Q.all(promises);
        })
        .then(function(){
            return Clotho.get([publicObjectId, privateObjectId, publicFunctionId, privateFunctionId]);
        })
        .then(function(result){
            objectCopies = result;
            return Clotho.logout();
        })
        .finally(function(){
            done();
        });
    });

    afterAll(function(done){
        Clotho.login("owner", "owner")
        .then(function(){
            return Clotho.destroy([publicObjectId,privateObjectId,publicFunctionId,privateFunctionId]);
        })
        .then(function(){
            return Clotho.logout();
        })
        .finally(function(){
            done();
        });
    });

    afterEach(function(done){
        Clotho.logout()
        .finally(function(){
            done();
        });
    });

    var testPublicInteractions = function(){
        testRead(publicObject, true);
        testRun(publicFunction, true);
    }
    
    var testRead = function(idString, expectSuccess){
        it((expectSuccess? "can" : "cannot") + " read that object", function(done){
            var id = resolveId(idString);
            Clotho.get(id)
            .then(function(result){
                if(expectSuccess) 
                    expect(result).not.toBeUndefined();
                else
                    expect(result).toBeUndefined();
            })
            .finally(function(){
                done();
            });
        })
    }
    
    var testRun = function(idString, expectSuccess){
        it((expectSuccess? "can" : "cannot") + " run that object", function(done){
            var id = resolveId(idString);
            Clotho.run({"module":id, function:"function"})
            .then(function(result){
                if(expectSuccess) 
                    expect(result).toEqual("function ran!");
                else
                    expect(result).toBeUndefined();
            })
            .finally(function(){
                done();
            });
        })
    }

    var testEdit = function(idString, expectSuccess){
        it((expectSuccess? "can" : "cannot") + " edit that object", function(done){
            var id = resolveId(idString);
            Clotho.set({"id":id, "newfield":"object edited"})
            .then(function(){
                return Clotho.logout();
            })
            .then(function(){
                return Clotho.login("owner", "owner");
            })
            .then(function(result){
                return Clotho.get(id);
            })
            .then(function(result){
                if (expectSuccess)
                    expect(result.newfield).toEqual("object edited");
                else 
                    expect(result.newfield).toBeUndefined();
            })
            .finally(function(){
                return Clotho.logout();
            })
            .finally(function(){
                return Clotho.login("owner", "owner");
            })
            .finally(function(){
                return Clotho.destroy(id);
            })
            .finally(function(){
                var copy;
                for (var i=0; i<objectCopies.length; i++){
                   if (objectCopies[i].id == id)
                        copy = objectCopies[i];
                }

                return Clotho.create({"id":copy.id, name:copy.name, schema:copy.schema});
            })
            .finally(function(){
                return Clotho.logout();
            })
            .finally(function(){
                done();
            });
        });
    }

    var testCreate = function(expectSuccess){
        it((expectSuccess? "can" : "cannot") + " create an object", function(done){
            Clotho.create({})
            .then(function(result){
                if(expectSuccess)
                    expect(result).not.toBeUndefined();
                else
                    expect(result).toBeUndefined();
            })
            .finally(function(){
                done();
            });
        });
    }

    var testGrant = function(idString, expectSuccess){
        it((expectSuccess? "can" : "cannot") + " grant privileges on that object", function(done){
            var id = resolveId(idString);
            Clotho.grant("none", id, ["write"], [])
            .then(function(result){
                return Clotho.logout();
            })
            .then(function(){
                return Clotho.login("owner", "owner");
            })
            .then(function(){
                return Clotho.get(id);
            })
            .then(function(result){
                if(expectSuccess)
                    expect(result.$$permissions.user.none).toContain("edit");
                else 
                    expect(result.$$permissions.user.none).toBeUndefined();
            })
            .finally(function(){
                return Clotho.logout();
            })
            .finally(function(){
                return Clotho.login("owner","owner");
            })
            .finally(function(result){
                //reset granted privileges
                Clotho.grant("none", id, [], ["write"]);
            })
            .finally(function(){
                done();
            });
        });
    }

    var testDelete = function(idString, expectSuccess){
        it((expectSuccess? "can" : "cannot") + " delete that object", function(done){
            var id = resolveId(idString);
            Clotho.destroy(id)
            .then(function(){
                return Clotho.logout();
            })
            .then(function(){
                return Clotho.login("owner", "owner");
            })
            .then(function(){
                return Clotho.get(id);
            })
            .then(function(result){
                if(expectSuccess)
                    expect(result).toBeUndefined();
                else
                    expect(result).not.toBeUndefined();
            })
            .finally(function(){
                var copy;
                for (var i=0; i<objectCopies.length; i++){
                   if (objectCopies[i].id == id)
                        copy = objectCopies[i];
                }

                return Clotho.create({"id":copy.id, name:copy.name, schema:copy.schema});
            })
            .finally(function(){
                done();
            });
        });
    }

    var loginAs = function(uname, pass){
        return function(done){
            Clotho.login(uname, pass)
            .finally(function(){
                done();
            });
        }
    }

    describe("when unauthenticated, a user", function(){
        //no login
        testCreate(false);

        describe("interacting with a public object", function(){
            testPublicInteractions();
            testEdit(publicObject, false);
            testGrant(publicObject, false);
            testDelete(publicObject, false);
        });
        
        describe("interacting with a private object", function(){
            testRead(privateObject, false);
            testRun(privateFunction, false);
            testEdit(privateObject, false);
            testGrant(privateObject, false);
            testDelete(privateObject, false);
        });
    });

    describe("an authenticated user", function(){
        describe("with no privileges", function(){ 
            beforeEach(loginAs("none", "none"));

            testCreate(true);

            describe("on a public object", function(){
                testPublicInteractions();
                testEdit(publicObject, false);
                testGrant(publicObject, false);
                testDelete(publicObject, false);
            });

            describe("on a private object", function(){
                testRead(privateObject, false);
                testRun(privateFunction, false);
                testEdit(privateObject, false);
                testGrant(privateObject, false);
                testDelete(privateObject, false);
            });
        });

        describe("with run privileges", function(){ 
            beforeEach(loginAs("run", "run"));

            testCreate(true);

            describe("on a public object", function(){
                testPublicInteractions();
                testEdit(publicObject, false);
                testGrant(publicObject, false);
                testDelete(publicObject, false);
            });

            describe("on a private object", function(){
                testRead(privateObject, false);
                testRun(privateFunction, true);
                testEdit(privateObject, false);
                testGrant(privateObject, false);
                testDelete(privateObject, false);
            });

        });

        describe("with read privileges", function(){ 
            beforeEach(loginAs("read", "read"));

            testCreate(true);

            describe("on a public object", function(){
                testPublicInteractions();
                testEdit(publicObject, false);
                testGrant(publicObject, false);
                testDelete(publicObject, false);
            });

            describe("on a private object", function(){
                testRead(privateObject, true);
                testRun(privateFunction, true);
                testEdit(privateObject, false);
                testGrant(privateObject, false);
                testDelete(privateObject, false);
            });
        });

        describe("with edit privileges", function(){
            beforeEach(loginAs("write", "write"));

            testCreate(true);

            describe("on a public object", function(){
                testPublicInteractions();
                testEdit(publicObject, true);
                testGrant(publicObject, false);
                testDelete(publicObject, false);
            });

            describe("on a private object", function(){
                testRead(privateObject, true);
                testRun(privateFunction, true);
                testEdit(privateObject, true);
                testGrant(privateObject, false);
                testDelete(privateObject, false);
            });
        });

        describe("with owner privileges", function(){
            beforeEach(loginAs("owner", "owner"));

            testCreate(true);

            describe("on a public object", function(){
                testPublicInteractions();
                testEdit(publicObject, true);
                testGrant(publicObject, true);
                testDelete(publicObject, true);
            });

            describe("on a private object", function(){
                testRead(privateObject, true);
                testRun(privateFunction, true);
                testEdit(privateObject, true);
                testGrant(privateObject, true);
                testDelete(privateObject, true);
            });
        });
    });
/*
    var passedFunction = function(done){
        expect(true).toBe(true);
        done();
    }

    var generatedFunction = function(i){
        return function(done){
            expect(i).toBe(1);
            done();
        }
    }
    describe("programmatically generating asynchronous tests", function(){
        it("passing in an async function", passedFunction);

        it("inline function", function(done){
            passedFunction(done);
        });
        
        it("normal", function(done){
            done();
        });

        it("passing in a generated async function", generatedFunction(1));
    });
    */
});
