/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.clothocad.core.security;

import com.fasterxml.jackson.annotation.JsonCreator;
import com.fasterxml.jackson.annotation.JsonProperty;

/**
 *
 * @author spaige
 */
public class DummyAccount extends UnauthenticableAccount {

    
    @JsonCreator
    public DummyAccount( @JsonProperty("_id") String username) {
        super(username);
    }

    @Override
    public void addGroup(String group) {
    }

    @Override
    public void addPermission(String permission) {
    }

    @Override
    public void removeGroup(String group) {
    }

    @Override
    public void removePermission(String permission) {
    }

}
